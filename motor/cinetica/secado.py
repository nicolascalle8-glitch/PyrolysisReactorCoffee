"""
secado.py — Evaporación de la humedad del lecho.

Recibe: agua (kg/m³ por nodo), T (K), dt (s).
Devuelve: (agua_nueva, evaporada, q) con q el calor promedio del paso en W/m³ (negativo: evaporar enfría).
Modelo de legacy [Supuesto]: primer orden de Arrhenius, k = A·exp(−Ea/RT) solo si T > 350 K, con
A = 5.0e5 1/s, Ea = 41.8 kJ/mol, ΔH = −2.26e6 J/kg. Integración exacta: evaporada = agua·(1 − e^(−k·dt)),
así nunca se evapora más agua de la que hay y la masa cierra para cualquier dt.
"""

import numpy as np

R_UNIVERSAL = 8.31446   # J/(mol·K)
A_EVAP = 5.0e5          # 1/s
EA_EVAP = 41.8e3        # J/mol
DH_EVAP = -2.26e6       # J/kg (signo: positivo = exotérmico)
T_UMBRAL = 350.0        # K


def paso_secado(agua, T, dt):
    T = np.maximum(T, 273.15)
    k = np.where(T > T_UMBRAL, A_EVAP * np.exp(-EA_EVAP / (R_UNIVERSAL * T)), 0.0)
    evaporada = agua * -np.expm1(-k * dt)
    return agua - evaporada, evaporada, evaporada / dt * DH_EVAP
