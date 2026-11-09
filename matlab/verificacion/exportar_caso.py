"""
exportar_caso.py — Caso térmico sin reacción para verificar el motor en MATLAB (prueba 7 del spec).

Escribe en la carpeta de salida:
  caso_conduccion.csv       un renglón por nodo: i, k, r, z, V, k_W_mK, rho_cp, T0, T_<t1>, T_<t2>, ...
                            (T de Python con Euler implícito en cada tiempo, en K)
  caso_conduccion_meta.csv  escalares: R_lecho, R_ext, dr_lecho, dr_acero, dz, L, Nr, Nl, Nz, h_ext,
                            T_horno, R_contacto_extra, T0, dt_python
Caso: malla base, lecho virgen a 25 °C con k y ρcp constantes (evaluados a 298.15 K), acero de base.toml,
horno fijo a 500 °C, h = h_c de base.toml, 30 min. Uso: python matlab/verificacion/exportar_caso.py
"""

import csv
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(AQUI))
sys.path.insert(0, RAIZ)

from motor import cargar_parametros  # noqa: E402
from motor.conduccion import conductancias, paso_conduccion  # noqa: E402
from motor.malla import crear_malla  # noqa: E402
from motor.propiedades import cp_lecho, k_lecho  # noqa: E402

T0, T_HORNO, T_FIN, DT, CADA = 298.15, 773.15, 1800.0, 0.5, 120.0


def _caso():
    p = cargar_parametros(os.path.join(RAIZ, "parametros", "base.toml"))
    m = crear_malla(p)
    k = np.full(m.forma, p["acero.k"])
    rc = np.full(m.forma, p["acero.rho"] * p["acero.cp"])
    k[:m.n_lecho] = k_lecho(T0, 0.0)
    rc[:m.n_lecho] = p["lecho.rho_aparente"] * cp_lecho(T0, 0.0)
    h = np.full(m.forma[1], p["horno.h_c"])
    Rc = np.full(m.forma[1], p["lecho.R_contacto_extra"])
    return p, m, k, rc, h, Rc


def _python(m, k, rc, h, Rc):
    T = np.full(m.forma, T0)
    tiempos, campos = [], []
    n = int(round(T_FIN / DT))
    cada = int(round(CADA / DT))
    for j in range(1, n + 1):
        T, _ = paso_conduccion(m, T, k, rc, h, T_HORNO, Rc, DT)
        if j % cada == 0:
            tiempos.append(j * DT); campos.append(T.copy())
    return tiempos, campos


def exportar(carpeta=AQUI):
    p, m, k, rc, h, Rc = _caso()
    tiempos, campos = _python(m, k, rc, h, Rc)
    os.makedirs(carpeta, exist_ok=True)
    nodos = os.path.join(carpeta, "caso_conduccion.csv")
    with open(nodos, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["i", "k", "r", "z", "V", "k_W_mK", "rho_cp", "T0"] + [f"T_{t:g}" for t in tiempos])
        for i in range(m.forma[0]):
            for kz in range(m.forma[1]):
                w.writerow([i, kz, f"{m.r[i]:.12g}", f"{m.z[kz]:.12g}", f"{m.V[i, kz]:.12g}", f"{k[i, kz]:.12g}",
                            f"{rc[i, kz]:.12g}", T0] + [f"{c[i, kz]:.10f}" for c in campos])
    meta = os.path.join(carpeta, "caso_conduccion_meta.csv")
    with open(meta, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["nombre", "valor"])
        for nombre, valor in (("R_lecho", m.R_lecho), ("R_ext", m.R_ext), ("dr_lecho", m.dr_lecho),
                              ("dr_acero", m.dr_acero), ("dz", m.dz), ("L", p["geometria.L"]),
                              ("Nr", m.forma[0]), ("Nl", m.n_lecho), ("Nz", m.forma[1]), ("h_ext", h[0]),
                              ("T_horno", T_HORNO), ("R_contacto_extra", Rc[0]), ("T0", T0), ("dt_python", DT)):
            w.writerow([nombre, f"{valor:.12g}"])
    return nodos, meta


def verificar_con_bdf(carpeta=AQUI):
    """Segundo camino en Python: C·dT/dt = −A·T + b integrado con BDF (scipy), igual que ode15s en MATLAB."""
    from scipy.integrate import solve_ivp
    _, m, k, rc, h, Rc = _caso()
    tiempos, campos = _python(m, k, rc, h, Rc)
    Nr, Nz = m.forma
    G_r, G_z = conductancias(m, k, Rc)
    G_ext = m.A_r[-1] * m.dz_k * h
    C = (rc * m.V).ravel()

    def f(_, y):
        T = y.reshape(Nr, Nz)
        q = np.zeros_like(T)
        flujo_r = G_r * (T[1:] - T[:-1]); q[:-1] += flujo_r; q[1:] -= flujo_r
        flujo_z = G_z * (T[:, 1:] - T[:, :-1]); q[:, :-1] += flujo_z; q[:, 1:] -= flujo_z
        q[-1] += G_ext * (T_HORNO - T[-1])
        return q.ravel() / C

    sol = solve_ivp(f, (0, T_FIN), np.full(Nr * Nz, T0), method="BDF", t_eval=tiempos, rtol=1e-8, atol=1e-6)
    err = max(np.max(np.abs(sol.y[:, j].reshape(Nr, Nz) - campos[j])) for j in range(len(tiempos)))
    return {"error_max_K": float(err), "dT0_K": T_HORNO - T0}


if __name__ == "__main__":
    for ruta in exportar():
        print("escrito:", ruta)
    print("segundo camino (BDF en Python):", verificar_con_bdf())
