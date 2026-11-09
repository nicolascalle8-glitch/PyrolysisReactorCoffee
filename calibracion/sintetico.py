"""
sintetico.py — Prueba del ajuste con datos sintéticos (subproyecto B, paso 4).

Genera las dos corridas con valores conocidos (VERDAD), les pone ruido gaussiano de 1 K, corre ajuste_termico.py
y compara lo recuperado con la verdad. Los datos NO van al repo: quedan en 01_Horno_Vertical/resultados/.

Corrida 1 (tubo vacío): el motor no tiene tubo vacío (decisión 10), así que se arma aquí con las piezas del motor:
  conducción 2D en el acero (paso_conduccion, con los nodos del lecho sin masa ni conducción), HornoFisico, PID, SSR y
  termopar de control. Programa del spec §4 sin salida manual: rampa 10 °C/min a 550 °C con PID, 20 min de meseta y
  enfriamiento libre con el SSR apagado. El ajuste usa un modelo 0D distinto al que genera los datos.
Corrida 2 (café): motor en modo medido con T_horno(t) prescrita (rampa, meseta, enfriamiento exponencial), con malla y
  pasos más finos que los del ajuste, para que el ajuste no use exactamente el mismo modelo numérico que generó los datos.

En la Terminal de macOS:
  .venv/bin/python calibracion/sintetico.py
"""

import csv
import json
import os
import sys
import time

import numpy as np

CARPETA = os.path.dirname(os.path.abspath(__file__))
if CARPETA not in sys.path:
    sys.path.insert(0, CARPETA)

import ajuste_termico as aj  # noqa: E402
from motor import cargar_parametros, simular  # noqa: E402
from motor.conduccion import paso_conduccion  # noqa: E402
from motor.control import PID, SSR, Programa  # noqa: E402
from motor.horno import HornoFisico  # noqa: E402
from motor.malla import crear_malla  # noqa: E402
from motor.termopar import Termopar  # noqa: E402

# Verdad del ensayo sintético: distinta del arranque de base.toml a propósito. Elegida para que P_max = 2500 W alcance
# (P_max/C_h = 12.5 K/min, UA_p·ΔT = 1575 W a 550 °C). No son propiedades del horno real.
VERDAD = {"horno.C_h": 12000.0, "horno.UA_p": 3.0, "horno.h_c": 35.0, "horno.emisividad": 0.6,
          "lecho.k_virgen_0": 0.15, "lecho.k_char_0": 0.11, "lecho.R_contacto_extra": 0.015}
RUIDO_K = 1.0
TAU_ENFRIAMIENTO = 4000.0   # s, solo para la T_horno prescrita de la corrida 2
NUMERICA_FINA = {"malla.Nr_lecho": 28, "malla.Nr_acero": 6, "malla.Nz": 5, "numerica.dt0": 1.0,
                 "numerica.dt_max": 2.0, "numerica.dT_max": 1.0, "numerica.dX_max": 0.02, "numerica.t_registro": 1.0}
COLUMNAS = ("t_s", "T_horno_C", "T_pared_C", "T_centro_C", "T_control_C", "T_set_C", "u_ssr_pct", "T_amb_C")
SALIDA = os.path.join(aj.RAIZ, "..", "..", "resultados", "calibracion_sintetica")


def _escribir(ruta, filas):
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(COLUMNAS)
        w.writerows(filas)


def _texto(v):
    return "" if v is None else f"{v:.4f}"


def generar_corrida1(ruta, semilla=1, ruido=RUIDO_K):
    """Tubo vacío con PID y dt = 1 s. Registra al inicio de cada paso el u que se aplica en ese paso."""
    p = cargar_parametros(aj.BASE, {**VERDAD, "horno.modo": "fisico", "programa.T_final": 550.0, "programa.meseta": 20.0})
    m = crear_malla(p)
    Nl = m.n_lecho
    k = np.full(m.forma, p["acero.k"])
    rc = np.full(m.forma, p["acero.rho"] * p["acero.cp"])
    k[:Nl], rc[:Nl] = 1e-6, 1e-3     # tubo vacío: el lecho no guarda calor ni conduce
    R_extra = np.zeros(m.forma[1])
    prog, ssr = Programa(p), SSR(p["control.modo_ssr"], p["control.t_ciclo"])
    pid = PID(p["control.Kp"], p["control.Ki"], p["control.Kd"], p["control.t_muestreo"])
    T_amb = p["horno.T_amb"] + aj.C_A_K
    horno = HornoFisico(p["horno.C_h"], p["horno.UA_p"], T_amb)
    tc = Termopar(p["termopares.tau"])
    h_c, eps, P_max = p["horno.h_c"], p["horno.emisividad"], p["horno.P_max"]
    mid = m.forma[1] // 2
    T = np.full(m.forma, prog.T_ini)
    T_h = T_tc = prog.T_ini
    rng = np.random.default_rng(semilla)
    filas, t, dt = [], 0.0, 1.0
    while t <= prog.t_total or T_h > 100.0 + aj.C_A_K:
        u = pid.calcular(t, prog.setpoint(t), T_tc) if t <= prog.t_total else 0.0
        Tw = float(T[-1, mid])
        filas.append([f"{t:.0f}", _texto(T_h - aj.C_A_K + rng.normal(0, ruido)), _texto(Tw - aj.C_A_K + rng.normal(0, ruido)),
                      _texto(None), _texto(T_tc - aj.C_A_K + rng.normal(0, ruido)), _texto(prog.setpoint(t) - aj.C_A_K),
                      _texto(100.0 * u), _texto(T_amb - aj.C_A_K)])
        P = ssr.energia(t, t + dt, u, P_max) / dt
        h_ext = aj.h_externo(T_h, T[-1], h_c, eps)
        T, Q = paso_conduccion(m, T, k, rc, h_ext, T_h, R_extra, dt)
        T_h = horno.paso(T_h, P, Q / dt, dt)
        T_tc = tc.paso(T_tc, float(T[-1, mid]), dt)
        t += dt
    _escribir(ruta, filas)
    return len(filas)


def perfil_horno_corrida2(T_ini=25.0, rampa=10.0, T_fin=500.0, meseta=30.0, T_amb=25.0, T_para=90.0):
    """T_horno(t) prescrita a 1 s: rampa, meseta y enfriamiento exponencial hasta T_para (°C)."""
    t_rampa = (T_fin - T_ini) / rampa * 60.0
    t_cool = t_rampa + meseta * 60.0
    t_total = t_cool + TAU_ENFRIAMIENTO * np.log((T_fin - T_amb) / (T_para - T_amb))
    t = np.arange(0.0, np.floor(t_total) + 1.0)
    T = np.where(t < t_rampa, T_ini + rampa / 60.0 * t, T_fin)
    enf = t >= t_cool
    T[enf] = T_amb + (T_fin - T_amb) * np.exp(-(t[enf] - t_cool) / TAU_ENFRIAMIENTO)
    return t, T


def generar_corrida2(ruta, ruta_limpia, semilla=2, ruido=RUIDO_K):
    """Café en modo medido con numérica fina y la verdad. Devuelve la masa residual final (%)."""
    t, T_h = perfil_horno_corrida2()
    _escribir(ruta_limpia, [[f"{a:.0f}", f"{b:.4f}", "", "", "", "", "", "25.0"] for a, b in zip(t, T_h)])
    fijos = {c: VERDAD[c] for c in VERDAD if not c.startswith("horno.C_") and c != "horno.UA_p"}
    p = cargar_parametros(aj.BASE, {**fijos, **NUMERICA_FINA, "horno.modo": "medido", "horno.archivo_medido": ruta_limpia})
    res = simular(p)
    h = res.historia
    ts = h["t"] * 60.0
    Tp, Tc = np.interp(t, ts, h["T_pared"]), np.interp(t, ts, h["T_centro"])
    rng = np.random.default_rng(semilla)
    filas = [[f"{a:.0f}", f"{Th + rng.normal(0, ruido):.4f}", f"{p_ + rng.normal(0, ruido):.4f}",
              f"{c + rng.normal(0, ruido):.4f}", "", "", "", "25.0"] for a, Th, p_, c in zip(t, T_h, Tp, Tc)]
    _escribir(ruta, filas)
    return float(h["masa_residual"][-1])


def tabla(informe):
    """Filas (parámetro, verdad, recuperado, sigma, error %)."""
    filas = []
    for corrida in ("corrida1", "corrida2"):
        for c, v in informe[corrida]["valores"].items():
            ver = VERDAD[c]
            filas.append((c, ver, v, informe[corrida]["sigma"].get(c), 100.0 * (v - ver) / ver))
    return filas


def main(salida=SALIDA):
    salida = os.path.abspath(salida)
    os.makedirs(salida, exist_ok=True)
    c1, c2, limpia = (os.path.join(salida, n) for n in ("corrida1.csv", "corrida2.csv", "corrida2_T_horno_limpia.csv"))
    t0 = time.perf_counter()
    n1 = generar_corrida1(c1)
    masa = generar_corrida2(c2, limpia)
    t1 = time.perf_counter()
    informe = aj.correr(c1, c2, salida)
    t2 = time.perf_counter()
    filas = tabla(informe)
    informe["verdad"] = VERDAD
    informe["errores_pct"] = {f[0]: f[4] for f in filas}
    informe["generacion_s"], informe["ajuste_s"] = t1 - t0, t2 - t1
    informe["masa_residual_final_pct"] = masa
    with open(os.path.join(salida, "calibracion.json"), "w", encoding="utf-8") as f:
        json.dump(informe, f, ensure_ascii=False, indent=2)
    print(f"corrida 1: {n1} registros. Generación {t1 - t0:.0f} s, ajuste {t2 - t1:.0f} s.")
    print(f"{'parámetro':28s}{'verdad':>12s}{'recuperado':>13s}{'sigma':>11s}{'error %':>10s}")
    for c, v, r, s, e in filas:
        print(f"{c:28s}{v:12.5g}{r:13.5g}{(s if s is not None else float('nan')):11.2g}{e:10.2f}")
    print("RMS corrida 1:", informe["corrida1"]["rms"])
    print("RMS corrida 2:", informe["corrida2"]["rms"])
    print("Correlación R_contacto–k:", informe["correlacion_R_k"])
    print("Criterio de cierre:", informe["criterio_cierre"])
    print("Resultado en", salida)
    return informe


if __name__ == "__main__":
    main()
