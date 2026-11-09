"""
control.py — Programa de temperatura, PID discreto y SSR.

Programa: rampa lineal desde programa.T_inicial a programa.rampa hasta programa.T_final y luego meseta.
PID: u = Kp·e + Ki·∫e dt + Kd·de/dt, muestreado cada t_muestreo (entre muestras mantiene u), salida en
[0, 1]. Anti-windup por integración condicional: si u satura y el error empuja hacia la saturación, el
integrador no crece. Las unidades del MC-5438 (banda P en %, I y D en s) se agregan con el manual [Por validar].
SSR: 'promedio' entrega u·P_max todo el tiempo; 'onoff' prende P_max durante u·t_ciclo al inicio de cada
ventana de t_ciclo y apaga el resto (ventana de tiempo).
"""

import math

C_A_K = 273.15


class Programa:
    def __init__(self, p):
        self.T_ini = p["programa.T_inicial"] + C_A_K
        self.T_fin = p["programa.T_final"] + C_A_K
        self.beta = p["programa.rampa"] / 60.0                     # K/s
        self.t_rampa = max(self.T_fin - self.T_ini, 0.0) / self.beta
        self.t_total = self.t_rampa + p["programa.meseta"] * 60.0

    def setpoint(self, t: float) -> float:
        if t >= self.t_rampa:
            return self.T_fin
        return self.T_ini + self.beta * t


class PID:
    def __init__(self, Kp, Ki, Kd, t_muestreo):
        self.Kp, self.Ki, self.Kd, self.t_muestreo = Kp, Ki, Kd, t_muestreo
        self.integral = 0.0
        self._t_ultimo = None
        self._e_ultimo = 0.0
        self.u = 0.0

    def calcular(self, t, setpoint, medida) -> float:
        if self._t_ultimo is not None and t < self._t_ultimo + self.t_muestreo - 1e-9:
            return self.u
        e = setpoint - medida
        h = 0.0 if self._t_ultimo is None else t - self._t_ultimo
        derivada = (e - self._e_ultimo) / h if h > 0 else 0.0
        integral = self.integral + self.Ki * e * h
        u = self.Kp * e + integral + self.Kd * derivada
        saturado_arriba = u > 1.0 and e > 0
        saturado_abajo = u < 0.0 and e < 0
        if not (saturado_arriba or saturado_abajo):
            self.integral = integral
        self.u = min(max(self.Kp * e + self.integral + self.Kd * derivada, 0.0), 1.0)
        self._t_ultimo, self._e_ultimo = t, e
        return self.u


class SSR:
    def __init__(self, modo: str, t_ciclo: float):
        if modo not in ("promedio", "onoff"):
            raise ValueError(f"modo SSR '{modo}' no válido: promedio u onoff")
        self.modo, self.t_ciclo = modo, t_ciclo

    def potencia(self, t, u, P_max) -> float:
        if self.modo == "promedio":
            return u * P_max
        return P_max if (t % self.t_ciclo) < u * self.t_ciclo else 0.0

    def _tiempo_prendido(self, t, u):
        n = math.floor(t / self.t_ciclo)
        return n * u * self.t_ciclo + min(t - n * self.t_ciclo, u * self.t_ciclo)

    def energia(self, t0, t1, u, P_max) -> float:
        """Energía exacta entregada entre t0 y t1 con u constante (J)."""
        if self.modo == "promedio":
            return u * P_max * (t1 - t0)
        return P_max * (self._tiempo_prendido(t1, u) - self._tiempo_prendido(t0, u))
