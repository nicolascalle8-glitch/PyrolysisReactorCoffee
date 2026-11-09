"""
horno.py — Temperatura del horno que rodea el tubo.

HornoIdeal (como legacy): T_h sigue al setpoint con retardo de primer orden τ (piedra refractaria).
  Exacto en el paso con T_set constante: T_h = T_set + (T_h − T_set)·e^(−dt/τ).
HornoFisico: C_h·dT_h/dt = P − UA_p·(T_h − T_amb) − Q_tubo (spec §5.2), con P y Q_tubo constantes en el paso.
  Es lineal en T_h: se integra exacto. T∞ = T_amb + (P − Q_tubo)/UA_p; T_h = T∞ + (T_h − T∞)·e^(−UA_p·dt/C_h).
HornoMedido: T_h(t) leída de un CSV del ensayo (spec de B §5) e interpolada linealmente; el tiempo cuenta desde el
  primer registro. Columnas t_s y T_horno_C obligatorias; T_set_C y u_ssr_pct opcionales (NaN si faltan).
C_h, UA_p y la emisividad del intercambio son [Por validar]: se calibran en el subproyecto B.
"""

import csv
import math

import numpy as np

from .parametros import ErrorParametro

CLAVE_ARCHIVO = "horno.archivo_medido"
C_A_K = 273.15


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


def _numero(fila, columna, n_fila, ruta, obligatorio):
    """Valor numérico de una celda; una celda opcional vacía o no finita da NaN."""
    texto = (fila.get(columna) or "").strip()
    if not texto:
        if obligatorio:
            raise ErrorParametro(f"{CLAVE_ARCHIVO}: fila {n_fila} de {ruta}: falta el dato de {columna}")
        return math.nan
    try:
        valor = float(texto)
    except ValueError:
        raise ErrorParametro(f"{CLAVE_ARCHIVO}: fila {n_fila} de {ruta}: {columna} = '{texto}' no es un número") from None
    if not math.isfinite(valor):
        if obligatorio:
            raise ErrorParametro(f"{CLAVE_ARCHIVO}: fila {n_fila} de {ruta}: {columna} = '{texto}' no es finito")
        return math.nan
    return valor


class HornoMedido:
    def __init__(self, ruta):
        try:
            with open(ruta, newline="", encoding="utf-8-sig") as f:
                lector = csv.DictReader(f)
                columnas = lector.fieldnames or []
                for c in ("t_s", "T_horno_C"):
                    if c not in columnas:
                        raise ErrorParametro(f"{CLAVE_ARCHIVO}: a {ruta} le falta la columna {c}")
                t, T_h, T_set, u = [], [], [], []
                for fila in lector:
                    n = lector.line_num
                    t.append(_numero(fila, "t_s", n, ruta, True))
                    T_h.append(_numero(fila, "T_horno_C", n, ruta, True))
                    T_set.append(_numero(fila, "T_set_C", n, ruta, False))
                    u.append(_numero(fila, "u_ssr_pct", n, ruta, False))
                    if len(t) > 1 and t[-1] <= t[-2]:
                        raise ErrorParametro(f"{CLAVE_ARCHIVO}: fila {n} de {ruta}: t_s debe ser estrictamente creciente")
        except (OSError, UnicodeDecodeError) as e:
            raise ErrorParametro(f"{CLAVE_ARCHIVO}: no se pudo leer '{ruta}': {e}") from None
        if len(t) < 2:
            raise ErrorParametro(f"{CLAVE_ARCHIVO}: {ruta} necesita al menos 2 registros")
        self._t = np.array(t) - t[0]
        self._T_h = np.array(T_h) + C_A_K
        self._T_set = np.array(T_set) + C_A_K
        self._u = np.array(u) / 100.0
        self.t_fin = float(self._t[-1])

    def temperatura(self, t):
        """T del horno (K) en t (s desde el primer registro)."""
        return float(np.interp(t, self._t, self._T_h))

    def setpoint(self, t):
        """Setpoint registrado (K); NaN si el CSV no lo trae."""
        return float(np.interp(t, self._t, self._T_set))

    def salida(self, t):
        """Salida del SSR registrada (fracción 0–1); NaN si el CSV no la trae."""
        return float(np.interp(t, self._t, self._u))
