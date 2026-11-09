"""
horno.py — Temperatura del horno que rodea el tubo.

HornoIdeal (como legacy): T_h sigue al setpoint con retardo de primer orden τ (piedra refractaria).
  Exacto en el paso con T_set constante: T_h = T_set + (T_h − T_set)·e^(−dt/τ).
HornoFisico: C_h·dT_h/dt = P − UA_p·(T_h − T_amb) − Q_tubo (spec §5.2), con P y Q_tubo constantes en el paso.
  Es lineal en T_h: se integra exacto. T∞ = T_amb + (P − Q_tubo)/UA_p; T_h = T∞ + (T_h − T∞)·e^(−UA_p·dt/C_h).
C_h, UA_p y la emisividad del intercambio son [Por validar]: se calibran en el subproyecto B.
"""

import math


class HornoIdeal:
    def __init__(self, tau):
        self.tau = tau

    def paso(self, T_h, T_set, dt):
        return T_set + (T_h - T_set) * math.exp(-dt / self.tau)


class HornoFisico:
    def __init__(self, C_h, UA_p, T_amb):
        self.C_h, self.UA_p, self.T_amb = C_h, UA_p, T_amb

    def paso(self, T_h, P, Q_tubo, dt):
        neto = P - Q_tubo
        if self.UA_p == 0.0:
            return T_h + neto * dt / self.C_h
        T_inf = self.T_amb + neto / self.UA_p
        return T_inf + (T_h - T_inf) * math.exp(-self.UA_p * dt / self.C_h)
