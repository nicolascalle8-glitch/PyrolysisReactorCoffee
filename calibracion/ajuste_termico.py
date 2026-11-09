"""
ajuste_termico.py — Calibración térmica del horno y del lecho (subproyecto B, spec §4 y §5).

Corrida 1 (tubo vacío), decisión 10: modelo 0D del tubo, sin tocar el motor.
  C_tubo·dT_pared/dt = h_ext·A·(T_horno − T_pared),  h_ext = h_c + ε·σ·(T_h² + T_p²)(T_h + T_p)  (el mismo del motor)
  → h_c y ε con T_horno medido como entrada y T_pared como objetivo.
  C_h·dT_horno/dt = u·P_max − UA_p·(T_horno − T_amb) − Q_tubo  (el mismo balance de motor/horno.py)
  → C_h y UA_p con u_ssr_pct como entrada y T_horno como objetivo.
Corrida 2 (café), decisión 11: motor en horno.modo = "medido", con h_c y ε fijos de la corrida 1.
  → k_virgen_0 (25–200 °C), k_char_0 (enfriamiento final) y R_contacto_extra (débil: se reporta su correlación con k).
  La pirólisis (200 °C → fin de la meseta) no se ajusta; solo se compara (criterio, decisión 9).

Uso del motor: solo la interfaz pública (cargar_parametros, simular, simular_lote).
C_tubo y A salen de la geometría de base.toml: hay que poner la medida real antes de ajustar (spec §3.1).
La incertidumbre (sigma, correlación) sale del jacobiano: supone ruido independiente, así que con error de
modelo real es optimista. En la corrida 2, "sigma" suma además la incertidumbre de h_c y ε (propagación lineal);
"sigma_ajuste" y las correlaciones son solo del ajuste.

En la Terminal de macOS (uso real):
  .venv/bin/python calibracion/ajuste_termico.py --corrida1 c1.csv --corrida2 c2.csv --salida carpeta
"""

import argparse
import csv
import json
import os
import sys

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.optimize import least_squares

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if RAIZ not in sys.path:
    sys.path.insert(0, RAIZ)

from motor import cargar_parametros, simular, simular_lote  # noqa: E402

BASE = os.path.join(RAIZ, "parametros", "base.toml")
SIGMA = 5.670374419e-8   # W/(m²·K⁴), Stefan–Boltzmann (CODATA 2018), igual que el motor
C_A_K = 273.15
T_PIROLISIS = 200.0      # °C, spec §4: desde aquí no se ajusta
MARGEN_MESETA = 5.0      # K bajo T_máx del horno: de aquí en adelante empieza el enfriamiento
VENTANA_SUAVIZADO = 31   # muestras, solo para ubicar los tramos (no toca los datos del ajuste)
RMS_PARED, RMS_HORNO, RMS_CENTRO_FUERA, RMS_CENTRO_PIROLISIS = 5.0, 5.0, 10.0, 20.0  # K, decisión 9
CORRELACION_ALTA = 0.9   # |corr(R_contacto, k)| mayor que esto: dejar R_contacto fijo (spec §4)

# Numérica del motor durante el ajuste. Con caras z adiabáticas y T_horno uniforme la solución no depende de z,
# así que 3 nodos axiales dan lo mismo que 20 y corren ~7 veces más rápido. El paso de tiempo sí importa: con los
# valores de base.toml (dt_max = 20 s) el ajuste sintético sesgaba k y R_contacto hasta −34 %.
NUMERICA_AJUSTE = {"malla.Nz": 3, "numerica.t_registro": 10.0, "numerica.dt0": 1.0, "numerica.dt_max": 2.0,
                   "numerica.dT_max": 1.0, "numerica.dX_max": 0.02}
TOL_INICIO = 5.0         # K, diferencia máxima entre T_pared o T_centro y T_horno al arrancar
ANCHO_MESETA = 2.0       # K, ancho de las clases para ubicar la T de meseta
PASOS_FIJOS = {"horno.h_c": 1.0, "horno.emisividad": 0.02}   # pasos de la derivada respecto a los parámetros fijos


# ------------------------------------------------------------------ datos
def leer_ensayo(ruta, requeridas):
    """Lee el CSV del ensayo (spec §3.3). Devuelve dict de arreglos con t relativo al primer registro."""
    with open(ruta, newline="", encoding="utf-8") as f:
        filas = list(csv.DictReader(f))
    if len(filas) < 2:
        raise ValueError(f"{ruta}: tiene menos de 2 registros")
    columnas = list(filas[0])
    faltan = [c for c in requeridas if c not in columnas]
    if faltan:
        raise ValueError(f"{ruta}: faltan las columnas {faltan}")
    datos = {}
    for c in columnas:
        valores = np.empty(len(filas))
        for i, fila in enumerate(filas):
            texto = (fila[c] or "").strip()
            valores[i] = float(texto) if texto else np.nan
            if c in requeridas and not np.isfinite(valores[i]):
                raise ValueError(f"{ruta}: la columna {c} no tiene dato en la fila {i + 2}")
        datos[c] = valores
    datos["t_s"] = datos["t_s"] - datos["t_s"][0]
    if np.any(np.diff(datos["t_s"]) <= 0):
        raise ValueError(f"{ruta}: t_s debe ser estrictamente creciente")
    return datos


def _inicio_enfriamiento(ens, horno):
    """Último registro antes del enfriamiento libre: el último con SSR encendido; sin u_ssr_pct, el último cerca de la meseta.

    La T de meseta es la clase más poblada de la mitad alta del registro, no el máximo: así un sobrepaso no la mueve.
    """
    u = ens.get("u_ssr_pct")
    if u is not None and np.any(u > 0.0):
        return int(np.nonzero(u > 0.0)[0][-1])
    arriba = horno[horno >= 0.5 * (horno.min() + horno.max())]
    cuentas, bordes = np.histogram(arriba, bins=np.arange(arriba.min(), arriba.max() + ANCHO_MESETA, ANCHO_MESETA))
    k = int(np.argmax(cuentas))
    return int(np.nonzero(horno >= 0.5 * (bordes[k] + bordes[k + 1]) - MARGEN_MESETA)[0][-1])


def verificar_inicio(ens):
    """El motor arranca todo a T_horno(0): rechaza corridas donde pared o centro empiezan a otra temperatura."""
    T_h0 = float(np.mean(ens["T_horno_C"][:5]))
    for c in ("T_pared_C", "T_centro_C"):
        T0 = float(np.mean(ens[c][:5]))
        if abs(T0 - T_h0) > TOL_INICIO:
            raise ValueError(f"{c} arranca en {T0:.1f} °C y T_horno en {T_h0:.1f} °C (diferencia mayor que {TOL_INICIO} K): "
                             "el modelo supone tubo y lecho a T_horno(0); repetir la corrida desde el equilibrio")


def tramos(ens):
    """Separa la corrida 2 (spec §4): virgen (T_pared ≤ 200 °C al subir), enfriamiento y pirólisis (lo demás)."""
    pared = uniform_filter1d(ens["T_pared_C"], VENTANA_SUAVIZADO, mode="nearest")
    horno = uniform_filter1d(ens["T_horno_C"], VENTANA_SUAVIZADO, mode="nearest")
    n = pared.size
    i_max = int(np.argmax(horno))
    sobre = np.nonzero(pared[:i_max + 1] > T_PIROLISIS)[0]
    if sobre.size == 0:
        raise ValueError(f"La corrida nunca pasa de {T_PIROLISIS} °C en la pared")
    i_enf = _inicio_enfriamiento(ens, horno)
    virgen = np.arange(n) < sobre[0]
    enfriamiento = np.arange(n) > i_enf
    if virgen.sum() < 60 or enfriamiento.sum() < 60:
        raise ValueError("La corrida no tiene datos suficientes antes de 200 °C o en el enfriamiento")
    return {"virgen": virgen, "enfriamiento": enfriamiento, "pirolisis": ~(virgen | enfriamiento)}


# ------------------------------------------------------------------ utilidades del ajuste
def _rms(x):
    return float(np.sqrt(np.mean(np.square(x))))


def _incertidumbre(J, r, nombres):
    """Sigma y matriz de correlación con cov = s²·(JᵀJ)⁻¹; escala las columnas para no perder precisión."""
    n, m = J.shape
    escala = np.linalg.norm(J, axis=0)
    escala[escala == 0.0] = 1.0
    Js = J / escala
    cov_s = float(np.sum(r**2)) / max(n - m, 1) * np.linalg.pinv(Js.T @ Js)
    sig_s = np.sqrt(np.maximum(np.diag(cov_s), 0.0))
    divisor = np.outer(sig_s, sig_s)
    corr = np.divide(cov_s, divisor, out=np.zeros_like(cov_s), where=divisor > 0.0)
    return ({c: float(s) for c, s in zip(nombres, sig_s / escala)},
            {a: {b: float(corr[i, j]) for j, b in enumerate(nombres)} for i, a in enumerate(nombres)})


def sigma_con_fijos(J, cov_ajuste, J_fijos, cov_fijos):
    """Sigma del ajuste sumando la incertidumbre de parámetros fijos θ: dx/dθ = −J⁺·J_θ (propagación lineal)."""
    S = -np.linalg.lstsq(J, J_fijos, rcond=None)[0]
    return np.sqrt(np.maximum(np.diag(cov_ajuste + S @ cov_fijos @ S.T), 0.0))


def _covarianza(sigma, corr, nombres):
    s = np.array([sigma[c] for c in nombres])
    return np.outer(s, s) * np.array([[corr[a][b] for b in nombres] for a in nombres])


def _geometria_tubo(p):
    R, e, L = p["geometria.R_lecho"], p["geometria.e_acero"], p["geometria.L"]
    R_ext = R + e
    return 2.0 * np.pi * R_ext * L, p["acero.rho"] * p["acero.cp"] * np.pi * (R_ext**2 - R**2) * L


# ------------------------------------------------------------------ corrida 1: tubo vacío
def h_externo(T_h, T_p, h_c, eps):
    """Coeficiente horno–tubo del motor (T en K)."""
    return h_c + eps * SIGMA * (T_h**2 + T_p**2) * (T_h + T_p)


def tubo_0d(t, T_h, T_ini, h_c, eps, A, C_tubo):
    """Integra el tubo 0D sobre la malla de los datos. Euler implícito con h_ext del paso anterior (como el motor)."""
    T = np.empty(t.size)
    T[0] = T_p = T_ini
    th, dts = T_h.tolist(), np.diff(t).tolist()
    factor = eps * SIGMA
    for n, dt in enumerate(dts):
        Th = th[n + 1]
        a = (h_c + factor * (Th * Th + T_p * T_p) * (Th + T_p)) * A * dt / C_tubo
        T_p = (T_p + a * Th) / (1.0 + a)
        T[n + 1] = T_p
    return T


def horno_0d(t, T_h0, P, UA_p, C_h, T_amb, Q_tubo):
    """Balance del horno, exacto en cada paso con P y Q_tubo constantes (mismo método de motor/horno.py)."""
    T = np.empty(t.size)
    T[0] = T_h = T_h0
    for n, dt in enumerate(np.diff(t).tolist()):
        neto = P[n] - Q_tubo[n]
        if UA_p == 0.0:
            T_h = T_h + neto * dt / C_h
        else:
            T_inf = T_amb[n] + neto / UA_p
            T_h = T_inf + (T_h - T_inf) * np.exp(-UA_p * dt / C_h)
        T[n + 1] = T_h
    return T


def ajustar_corrida1(ens, p):
    """Devuelve dict con valores, sigma, correlación, RMS y h_c/ε. ens: leer_ensayo de la corrida 1."""
    A, C_tubo = _geometria_tubo(p)
    t = ens["t_s"]
    T_h = ens["T_horno_C"] + C_A_K
    T_p = ens["T_pared_C"] + C_A_K
    T_ini = float(np.mean(T_p[:5]))

    def r_tubo(x):
        return tubo_0d(t, T_h, T_ini, x[0], x[1], A, C_tubo) - T_p

    tubo = least_squares(r_tubo, [p["horno.h_c"], 0.3], bounds=([1.0, 0.0], [500.0, 1.0]), x_scale=[50.0, 0.3])
    h_c, eps = (float(v) for v in tubo.x)
    sig_t, corr_t = _incertidumbre(tubo.jac, tubo.fun, ("horno.h_c", "horno.emisividad"))

    T_pm = tubo_0d(t, T_h, T_ini, h_c, eps, A, C_tubo)
    Q_tubo = h_externo(T_h, T_pm, h_c, eps) * A * (T_h - T_pm)
    P = ens["u_ssr_pct"] / 100.0 * p["horno.P_max"]
    T_amb = ens["T_amb_C"] + C_A_K

    def r_horno(x):
        return horno_0d(t, T_h[0], P, x[1], x[0], T_amb, Q_tubo) - T_h

    horno = least_squares(r_horno, [p["horno.C_h"], p["horno.UA_p"]], bounds=([100.0, 0.0], [1e7, 500.0]),
                          x_scale=[1e4, 5.0])
    C_h, UA_p = (float(v) for v in horno.x)
    sig_h, corr_h = _incertidumbre(horno.jac, horno.fun, ("horno.C_h", "horno.UA_p"))

    return {"valores": {"horno.C_h": C_h, "horno.UA_p": UA_p, "horno.h_c": h_c, "horno.emisividad": eps},
            "sigma": {**sig_h, **sig_t},
            "correlacion": {"tubo": corr_t, "horno": corr_h},
            "rms": {"T_pared": _rms(tubo.fun), "T_horno": _rms(horno.fun)},
            "epsilon_significativa": bool(eps > 2.0 * sig_t["horno.emisividad"]),
            "A_m2": A, "C_tubo_J_K": C_tubo}


# ------------------------------------------------------------------ corrida 2: café con el motor
def ajustar_corrida2(ruta_csv, ens, fijos, base=BASE, ajustar_R=True, procesos=None, incert_fijos=None):
    """Ajusta k_virgen_0, k_char_0 y (opcional) R_contacto_extra con el motor en modo medido.

    fijos: {"horno.h_c": …, "horno.emisividad": …} de la corrida 1. El simulador corre completo en cada evaluación.
    incert_fijos: {"sigma": …, "corr": …} de esos fijos; si viene, "sigma" incluye su efecto y "sigma_ajuste" no.
    """
    verificar_inicio(ens)
    if not os.path.isabs(ruta_csv):
        ruta_csv = os.path.abspath(ruta_csv)
    nombres = ["lecho.k_virgen_0", "lecho.k_char_0"] + (["lecho.R_contacto_extra"] if ajustar_R else [])
    base_cambios = {"horno.modo": "medido", "horno.archivo_medido": ruta_csv, **NUMERICA_AJUSTE, **fijos}
    p0 = cargar_parametros(base, base_cambios)
    lim = {"lecho.k_virgen_0": (0.02, 1.0), "lecho.k_char_0": (0.02, 1.0), "lecho.R_contacto_extra": (0.0, 0.1)}
    escala = {"lecho.k_virgen_0": 0.1, "lecho.k_char_0": 0.1, "lecho.R_contacto_extra": 0.01}
    x0 = np.array([p0[c] for c in nombres])
    inf = np.array([lim[c][0] for c in nombres])
    sup = np.array([lim[c][1] for c in nombres])
    x_escala = np.array([escala[c] for c in nombres])

    t = ens["t_s"]
    bloques = tramos(ens)
    ajuste = bloques["virgen"] | bloques["enfriamiento"]

    def caso(x, extra=None):
        return cargar_parametros(base, {**base_cambios, **dict(zip(nombres, map(float, x))), **(extra or {})})

    def perfiles(res):
        h = res.historia
        ts = h["t"] * 60.0
        return np.interp(t, ts, h["T_pared"]), np.interp(t, ts, h["T_centro"])

    def residuo(res):
        Tp, Tc = perfiles(res)
        return np.concatenate([(Tp - ens["T_pared_C"])[ajuste], (Tc - ens["T_centro_C"])[ajuste]])

    memo = {}

    def f(x):
        clave = tuple(np.round(x, 12))
        if clave not in memo:
            memo[clave] = residuo(simular(caso(x)))
        return memo[clave]

    def jac(x):
        h = np.maximum(0.03 * np.abs(x), 0.03 * x_escala)
        h = np.where(x + h > sup, -h, h)
        casos = [caso(x + h[i] * np.eye(x.size)[i]) for i in range(x.size)]
        n_proc = procesos or min(x.size, os.cpu_count() or 1)
        resultados = simular_lote(casos, procesos=n_proc)
        for r in resultados:
            if isinstance(r, Exception):
                raise r
        r0 = f(x)
        return np.column_stack([(residuo(r) - r0) / h[i] for i, r in enumerate(resultados)])

    sol = least_squares(f, x0, jac=jac, bounds=(inf, sup), x_scale=x_escala, xtol=1e-3, ftol=1e-4, max_nfev=20)
    sigma, corr = _incertidumbre(sol.jac, sol.fun, nombres)
    sigma_ajuste = sigma
    if incert_fijos:
        n_fijos = list(fijos)
        pasos = [PASOS_FIJOS[c] if fijos[c] + PASOS_FIJOS[c] <= 1.0 or c != "horno.emisividad" else -PASOS_FIJOS[c]
                 for c in n_fijos]
        lote = simular_lote([caso(sol.x, {c: fijos[c] + h}) for c, h in zip(n_fijos, pasos)],
                            procesos=procesos or min(len(n_fijos), os.cpu_count() or 1))
        for r in lote:
            if isinstance(r, Exception):
                raise r
        r0 = f(sol.x)
        J_fijos = np.column_stack([(residuo(r) - r0) / h for r, h in zip(lote, pasos)])
        s_total = sigma_con_fijos(sol.jac, _covarianza(sigma, corr, nombres), J_fijos,
                                  _covarianza(incert_fijos["sigma"], incert_fijos["corr"], n_fijos))
        sigma = {c: float(s) for c, s in zip(nombres, s_total)}
    res_final = simular(caso(sol.x))
    Tp, Tc = perfiles(res_final)
    dp, dc = Tp - ens["T_pared_C"], Tc - ens["T_centro_C"]
    fuera = bloques["virgen"] | bloques["enfriamiento"]
    return {"valores": {c: float(v) for c, v in zip(nombres, sol.x)},
            "sigma": sigma, "sigma_ajuste": sigma_ajuste, "correlacion": corr,
            "rms": {"T_pared": _rms(dp), "T_centro_fuera_pirolisis": _rms(dc[fuera]),
                    "T_centro_pirolisis": _rms(dc[bloques["pirolisis"]])},
            "tramos_s": {n: [float(t[m][0]), float(t[m][-1])] for n, m in bloques.items()},
            "evaluaciones": int(sol.nfev),
            "avisos": list(res_final.avisos)}


# ------------------------------------------------------------------ cierre e informe
def criterio_cierre(c1, c2):
    """Decisión 9: RMS T_pared ≤ 5 K (dos corridas); T_horno ≤ 5 K (corrida 1); T_centro ≤ 10 K fuera de la pirólisis y ≤ 20 K en ella."""
    return {"T_pared_corrida1": c1["rms"]["T_pared"] <= RMS_PARED,
            "T_horno_corrida1": c1["rms"]["T_horno"] <= RMS_HORNO,
            "T_pared_corrida2": c2["rms"]["T_pared"] <= RMS_PARED,
            "T_centro_fuera_pirolisis": c2["rms"]["T_centro_fuera_pirolisis"] <= RMS_CENTRO_FUERA,
            "T_centro_pirolisis": c2["rms"]["T_centro_pirolisis"] <= RMS_CENTRO_PIROLISIS}


def correlacion_R_k(c2):
    """Correlación de R_contacto con cada k, y si conviene dejar R fijo (spec §4)."""
    corr = c2["correlacion"].get("lecho.R_contacto_extra")
    if corr is None:
        return None
    r_v, r_c = corr["lecho.k_virgen_0"], corr["lecho.k_char_0"]
    return {"k_virgen_0": r_v, "k_char_0": r_c, "dejar_R_fijo": bool(max(abs(r_v), abs(r_c)) > CORRELACION_ALTA)}


def correr(ruta_c1, ruta_c2, salida, base=BASE, ajustar_R=True):
    p = cargar_parametros(base)
    ens1 = leer_ensayo(ruta_c1, ("t_s", "T_horno_C", "T_pared_C", "u_ssr_pct", "T_amb_C"))
    c1 = ajustar_corrida1(ens1, p)
    ens2 = leer_ensayo(ruta_c2, ("t_s", "T_horno_C", "T_pared_C", "T_centro_C"))
    fijos = {c: c1["valores"][c] for c in ("horno.h_c", "horno.emisividad")}
    incert = {"sigma": {c: c1["sigma"][c] for c in fijos}, "corr": c1["correlacion"]["tubo"]}
    c2 = ajustar_corrida2(ruta_c2, ens2, fijos, base=base, ajustar_R=ajustar_R, incert_fijos=incert)
    informe = {"corrida1": c1, "corrida2": c2, "correlacion_R_k": correlacion_R_k(c2),
               "criterio_cierre": criterio_cierre(c1, c2)}
    os.makedirs(salida, exist_ok=True)
    with open(os.path.join(salida, "calibracion.json"), "w", encoding="utf-8") as f:
        json.dump(informe, f, ensure_ascii=False, indent=2)
    return informe


def main():
    a = argparse.ArgumentParser(description="Calibración térmica del horno y del lecho (subproyecto B)")
    a.add_argument("--corrida1", required=True, help="CSV de la corrida en tubo vacío")
    a.add_argument("--corrida2", required=True, help="CSV de la corrida con café")
    a.add_argument("--salida", required=True, help="carpeta donde queda calibracion.json")
    a.add_argument("--base", default=BASE, help="TOML de parámetros (con la geometría real)")
    a.add_argument("--R-fijo", action="store_true", help="deja R_contacto_extra en el valor del TOML")
    args = a.parse_args()
    informe = correr(args.corrida1, args.corrida2, args.salida, base=args.base, ajustar_R=not args.R_fijo)
    print(json.dumps(informe, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
