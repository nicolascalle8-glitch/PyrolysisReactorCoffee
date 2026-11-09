"""
tga.py — Muestra 0D térmicamente delgada (como un TGA): T impuesta, solo cambia la química.

simular_0d(regimen, rampa, T_final, meseta, humedad, cenizas=2.31, composicion=None) -> dict en % de la carga
(biochar, bio_oil, agua, gas, sin_reaccionar, humedad_restante).
Rampa en °C/min desde 25 °C, T_final en °C, meseta en min, humedad en % bh, cenizas en % bs.
Sirve para validar la cinética contra datos de literatura sin la transferencia de calor del reactor y,
en B, para ajustar la cinética contra un TGA. Mismo algoritmo que legacy/python/comparar_cineticas.simular_rac:
paso dt = min(1 s, 0.5 K / β), T al inicio del paso, integración exacta de la química y del secado.
"""

import numpy as np

from .rac import CineticaRAC
from .secado import paso_secado

T_AMB = 298.15  # K


def simular_0d(regimen, rampa, T_final, meseta, humedad, cenizas=2.31, composicion=None, acetilos_hce=0.1):
    beta = rampa / 60.0                       # K/s
    T_f = T_final + 273.15
    t_total = (T_f - T_AMB) / beta + meseta * 60.0
    dt = min(1.0, 0.5 / beta)
    kin = CineticaRAC(regimen, acetilos_hce=acetilos_hce)
    w, c = humedad / 100.0, cenizas / 100.0
    seco = 1.0 - w
    m = kin.estado_inicial((1.0 - c) * seco, composicion=composicion)   # base: 1 kg de café húmedo
    agua = np.array(w)
    acum = {"gas": 0.0, "agua_rxn": 0.0, "bio_oil": 0.0}
    vapor = 0.0
    t = 0.0
    while t < t_total:
        h = min(dt, t_total - t)
        T = min(T_AMB + beta * t, T_f)
        _, lib = kin.paso(m, T, h)
        for f in acum:
            acum[f] += float(lib[f])
        agua, ev, _ = paso_secado(agua, np.array(T), h)
        vapor += float(ev)
        t += h
    sin_reaccionar = float(kin.masa_original(m))
    r = {"biochar": float(kin.masa_solida(m)) - sin_reaccionar + c * seco, "bio_oil": acum["bio_oil"],
         "agua": acum["agua_rxn"] + vapor, "gas": acum["gas"], "sin_reaccionar": sin_reaccionar,
         "humedad_restante": float(agua)}
    return {k: 100.0 * v for k, v in r.items()}
