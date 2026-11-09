"""
simulador.py — Bucle de tiempo del motor (spec §5). No tiene física adentro: llama a los módulos.

Un paso dt:
 1. Control: setpoint del programa; en modo físico el PID lee el termopar de control y el SSR entrega
    la energía del paso.
 2. Horno: modo ideal = retardo al setpoint (antes de la conducción, como legacy); modo físico = balance
    de energía del horno con el calor real que salió hacia el tubo en este paso (después de la conducción),
    así horno y tubo intercambian exactamente la misma energía.
 3. Conducción implícita (lecho + acero) con k y ρcp del paso anterior; h_ext = h_c + h_rad linealizada.
 4. Química con Tⁿ⁺¹: RAC + TGL y secado, integración exacta; el calor de reacción ajusta T (ΔT = q·dt/ρcp).
 5. Control de dt (spec §5.5): si algún nodo cambia más de numerica.dT_max o la conversión más de
    numerica.dX_max, se repite el paso con dt/2; si el cambio es menor que la mitad del límite, dt crece ×1.5
    hasta numerica.dt_max. Con horno físico ningún paso cruza un instante de muestreo del PID
    (control.t_muestreo), y con SSR onoff dt ≤ t_ciclo. Piso DT_MIN para no quedar en bucle.
 6. Registro cada numerica.t_registro.
Balance de energía por paso: Σ ρcp·V·ΔT = Q_externo + Q_reacción (debe dar ~0 por construcción).
"""

import copy

import numpy as np

from .cinetica import CineticaRAC, paso_secado
from .conduccion import paso_conduccion
from .control import PID, SSR, Programa
from .horno import HornoFisico, HornoIdeal
from .malla import crear_malla
from .propiedades import coeficientes, cp_lecho, k_lecho, rho_lecho
from .resultados import Resultado
from .termopar import Termopar, nodo_mas_cercano

SIGMA = 5.670374419e-8   # W/(m²·K⁴), Stefan–Boltzmann (CODATA 2018)
C_A_K = 273.15
N_MAPAS = 12
DT_MIN = 1e-2            # s, piso del dt adaptativo (si se llega aquí, el paso se acepta y se cuenta)
CRECE, REDUCE = 1.5, 0.5


class ErrorSimulacion(RuntimeError):
    """La simulación no se puede seguir (T no finita o divergente)."""


class _Estado:
    """Todo lo que cambia en el tiempo, para poder repetir un paso (dt adaptativo)."""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def copia(self):
        nuevo = {}
        for c, v in self.__dict__.items():
            if isinstance(v, np.ndarray):
                nuevo[c] = v.copy()
            elif isinstance(v, dict):
                nuevo[c] = dict(v)
            else:
                nuevo[c] = v
        return _Estado(**nuevo)


class Simulacion:
    def __init__(self, p):
        self.p = p
        self.m = m = crear_malla(p)
        self.prog = Programa(p)
        self.kin = CineticaRAC(p["cinetica.regimen"], acetilos_hce=p["cinetica.acetilos_hce"])
        Nl, Nz = m.n_lecho, m.forma[1]
        self.Nl = Nl

        w = p["lecho.humedad"] / 100.0
        rho = p["lecho.rho_aparente"]
        self.rho_ceniza = rho * p["lecho.cenizas"] / 100.0 * (1.0 - w)
        self.rho_agua0 = rho * w
        self.rho_org0 = rho - self.rho_agua0 - self.rho_ceniza
        comp = {c: p[f"lecho.{c}"] / 100.0 for c in ("celulosa", "hemicelulosa", "lignina", "lipidos")}
        comp["cenizas"] = p["lecho.cenizas"] / 100.0
        especies = self.kin.estado_inicial(self.rho_org0, composicion=comp, forma=(Nl, Nz))

        self.V_l = m.V[:Nl]
        self.M0 = float(np.sum((self.rho_org0 + self.rho_agua0 + self.rho_ceniza) * self.V_l))
        self.rho_cp_acero = p["acero.rho"] * p["acero.cp"]
        self.k_acero = p["acero.k"]
        self.R_extra = np.full(Nz, p["lecho.R_contacto_extra"])
        self.h_c, self.eps = p["horno.h_c"], p["horno.emisividad"]
        self.fisico = p["horno.modo"] == "fisico"
        self.horno = (HornoFisico(p["horno.C_h"], p["horno.UA_p"], p["horno.T_amb"] + C_A_K)
                      if self.fisico else HornoIdeal(p["horno.tau_piedra"]))
        self.pid = PID(p["control.Kp"], p["control.Ki"], p["control.Kd"], p["control.t_muestreo"])
        self.ssr = SSR(p["control.modo_ssr"], p["control.t_ciclo"])
        self.P_max = p["horno.P_max"]
        self.tc = Termopar(p["termopares.tau"])
        self.nodo_tc = nodo_mas_cercano(m, p["termopares.control_r"], p["termopares.control_z"])
        self.mid = Nz // 2
        self.coef = coeficientes(p)
        self.secado = {c: p[f"secado.{c}"] for c in ("A", "Ea", "dH", "T_umbral")}

        T0 = self.prog.T_ini
        self.e = _Estado(t=0.0, T=np.full(m.forma, T0), especies=especies,
                         agua=np.full((Nl, Nz), self.rho_agua0), T_h=T0, T_tc=T0, T_set=T0, u=0.0,
                         lib={"bio_oil": 0.0, "gas": 0.0, "agua_rxn": 0.0, "evap": 0.0},
                         E_ext=0.0, E_rxn=0.0, err_E=0.0, E_abs=0.0)

    # ---------------------------------------------------------------- propiedades
    def _propiedades(self, e):
        Nl = self.Nl
        biomasa = self.kin.masa_original(e.especies)
        conv = np.clip(1.0 - (biomasa + e.agua) / (self.rho_org0 + self.rho_agua0), 0.0, 1.0)
        Tl = e.T[:Nl]
        k = np.full(e.T.shape, self.k_acero)
        rc = np.full(e.T.shape, self.rho_cp_acero)
        k[:Nl] = k_lecho(Tl, conv, self.coef)
        rho = rho_lecho(np.stack([self.kin.masa_solida(e.especies), np.full_like(Tl, self.rho_ceniza), e.agua]))
        rc[:Nl] = rho * cp_lecho(Tl, conv, self.coef)
        return k, rc, conv

    # ---------------------------------------------------------------- un paso
    def paso(self, e, dt):
        """Avanza el estado e (se modifica) un paso dt. Devuelve (ΔT máx, Δconversión máx)."""
        t = e.t
        T_set = self.prog.setpoint(t + 0.5 * dt)
        if self.fisico:
            e.u = self.pid.calcular(t, self.prog.setpoint(t), e.T_tc)
            P = self.ssr.energia(t, t + dt, e.u, self.P_max) / dt
        else:
            e.T_h = self.horno.paso(e.T_h, T_set, dt)
        k, rc, conv0 = self._propiedades(e)
        Tw = e.T[-1]
        h_ext = self.h_c + self.eps * SIGMA * (e.T_h**2 + Tw**2) * (e.T_h + Tw)
        T_ant = e.T
        T_new, Q = paso_conduccion(self.m, e.T, k, rc, h_ext, e.T_h, self.R_extra, dt)
        if not np.all(np.isfinite(T_new)):
            i, kz = np.argwhere(~np.isfinite(T_new))[0]
            raise ErrorSimulacion(f"La temperatura no es finita en t = {t + dt:.1f} s, nodo (r = {i}, z = {kz})")

        Nl = self.Nl
        Tq = np.maximum(T_new[:Nl], C_A_K)
        q_pir, lib = self.kin.paso(e.especies, Tq, dt)
        e.agua, evap, q_sec = paso_secado(e.agua, Tq, dt, **self.secado)
        q = q_pir + q_sec
        T_new[:Nl] += q * dt / rc[:Nl]
        Q_rxn = float(np.sum(q * dt * self.V_l))
        for f in ("bio_oil", "gas", "agua_rxn"):
            e.lib[f] += float(np.sum(lib[f] * self.V_l))
        e.lib["evap"] += float(np.sum(evap * self.V_l))

        dE = float(np.sum(rc * self.m.V * (T_new - T_ant)))
        e.err_E += dE - Q - Q_rxn
        e.E_ext += Q
        e.E_rxn += Q_rxn
        e.E_abs += abs(Q)
        if self.fisico:
            e.T_h = self.horno.paso(e.T_h, P, Q / dt, dt)
        if not np.all(np.isfinite(T_new)) or not np.isfinite(e.T_h):
            bad = np.argwhere(~np.isfinite(T_new))
            i, kz = bad[0] if bad.size else (-1, -1)
            raise ErrorSimulacion(f"La temperatura no es finita en t = {t + dt:.1f} s, nodo (r = {i}, z = {kz})")
        e.T = T_new
        e.T_tc = self.tc.paso(e.T_tc, e.T[self.nodo_tc], dt)
        e.T_set = T_set
        e.t = t + dt
        _, _, conv1 = self._propiedades(e)
        return float(np.max(np.abs(T_new - T_ant))), float(np.max(np.abs(conv1 - conv0)))

    # ---------------------------------------------------------------- balances y registro
    def _masas(self, e):
        V = self.V_l
        biomasa = float(np.sum(self.kin.masa_original(e.especies) * V))
        solido = float(np.sum(self.kin.masa_solida(e.especies) * V))
        ceniza = self.rho_ceniza * float(np.sum(V))
        agua = float(np.sum(e.agua * V))
        pct = lambda x: 100.0 * x / self.M0
        r = {"biochar": pct(solido - biomasa + ceniza), "bio_oil": pct(e.lib["bio_oil"]), "gas": pct(e.lib["gas"]),
             "agua": pct(e.lib["evap"] + e.lib["agua_rxn"]), "sin_reaccionar": pct(biomasa),
             "humedad_restante": pct(agua)}
        return r, pct(solido + ceniza + agua)

    def ejecutar(self, dt_fijo=None):
        p = self.p
        e = self.e
        t_fin = self.prog.t_total
        t_reg = p["numerica.t_registro"]
        hist = {c: [] for c in ("t", "T_horno", "T_set", "T_pared", "T_lecho_superficie", "T_centro",
                                "T_termopar_control", "masa_residual", "u")}
        mapas, err_masa, err_E = [], 0.0, 0.0
        t_mapas = list(np.linspace(0.0, t_fin, N_MAPAS))
        ultimo_reg = -np.inf
        pasos = 0

        def registrar():
            nonlocal err_masa, err_E
            r, residual = self._masas(e)
            err_masa = max(err_masa, abs(sum(r.values()) - 100.0))
            err_E = max(err_E, 100.0 * abs(e.err_E) / max(e.E_abs, 1e-30))
            T = e.T - C_A_K
            hist["t"].append(e.t / 60.0); hist["T_horno"].append(e.T_h - C_A_K); hist["T_set"].append(self.prog.setpoint(e.t) - C_A_K)
            hist["T_pared"].append(T[-1, self.mid]); hist["T_lecho_superficie"].append(T[self.Nl - 1, self.mid])
            hist["T_centro"].append(T[0, self.mid]); hist["T_termopar_control"].append(e.T_tc - C_A_K)
            hist["masa_residual"].append(residual); hist["u"].append(e.u)
            return r

        dt = dt_fijo or p["numerica.dt0"]
        dt_max = p["numerica.dt_max"]
        lim_T, lim_X = p["numerica.dT_max"], p["numerica.dX_max"]
        if self.ssr.modo == "onoff" and self.fisico:
            dt_max = min(dt_max, self.ssr.t_ciclo)
            dt = min(dt, dt_max)
        dT_max_aceptado, en_dt_min = 0.0, 0
        while e.t < t_fin - 1e-9:
            if e.t - ultimo_reg >= t_reg - 1e-9:
                registrar(); ultimo_reg = e.t
            while t_mapas and e.t >= t_mapas[0] - 1e-9:
                mapas.append((e.t / 60.0, e.T - C_A_K)); t_mapas.pop(0)
            if dt_fijo:
                self.paso(e, min(dt, t_fin - e.t))
                pasos += 1
                continue
            while True:
                h = min(dt, t_fin - e.t)
                if self.fisico:  # el PID muestrea cada t_muestreo: el paso no cruza la siguiente muestra
                    t_muestra = self.pid.siguiente_muestra()
                    if t_muestra is None or t_muestra <= e.t + 1e-9:   # muestrea al inicio de este paso
                        t_muestra = e.t + self.pid.t_muestreo
                    h = min(h, t_muestra - e.t)
                foto, pid_foto = e.copia(), copy.deepcopy(self.pid)
                dT, dX = self.paso(e, h)
                if (dT <= lim_T and dX <= lim_X) or dt <= DT_MIN:
                    break
                e, self.pid = foto, pid_foto
                self.e = e
                dt = max(dt * REDUCE, DT_MIN)
            if dT > lim_T or dX > lim_X:
                en_dt_min += 1
            else:
                dT_max_aceptado = max(dT_max_aceptado, dT)
            pasos += 1
            if dT < 0.5 * lim_T and dX < 0.5 * lim_X:
                dt = min(dt * CRECE, dt_max)
        r = registrar()
        if not mapas or mapas[-1][0] < e.t / 60.0 - 1e-9:
            mapas.append((e.t / 60.0, e.T - C_A_K))
        res = Resultado(rendimientos=r, historia={c: np.array(v) for c, v in hist.items()}, mapas=mapas,
                        balance_masa=err_masa, balance_energia=err_E, supuestos=p.supuestos(),
                        r_mm=self.m.r * 1e3, z_cm=self.m.z * 1e2, n_lecho=self.Nl, pasos=pasos)
        if not dt_fijo:
            res.dT_max_paso, res.pasos_en_dt_min = dT_max_aceptado, en_dt_min
        if not res.confiable:
            res.avisos.append(f"Balance fuera de tolerancia: masa {err_masa:.3g} %, energía {err_E:.3g} %")
        return res


def simular(p) -> Resultado:
    """Corre un caso completo con los parámetros p (ver spec §7)."""
    return Simulacion(p).ejecutar()
