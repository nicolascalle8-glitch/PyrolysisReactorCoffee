"""
termopar.py — Lectura de un termopar con retardo de primer orden.

Termopar(tau).paso(lectura, T_real, dt): exacto con T_real constante en el paso,
lectura = T_real + (lectura − T_real)·e^(−dt/τ). La posición (r, z) la escoge el simulador con
nodo_mas_cercano(malla, r, z). τ es [Supuesto] hasta medirlo.
"""

import math

import numpy as np


class Termopar:
    def __init__(self, tau):
        self.tau = tau

    def paso(self, lectura, T_real, dt):
        return T_real + (lectura - T_real) * math.exp(-dt / self.tau)


def nodo_mas_cercano(malla, r, z):
    return int(np.argmin(np.abs(malla.r - r))), int(np.argmin(np.abs(malla.z - z)))
