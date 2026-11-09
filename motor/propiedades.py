"""
propiedades.py — k, cp y ρ efectivas del lecho de café.

Recibe: T en K y la conversión (0 = café virgen, 1 = char) por nodo.
Devuelve: arreglos de la misma forma.
Fuente: correlaciones lineales de legacy/simulador_interactivo.py (k_bed, cp_bed, rho_bed), [Supuesto, sin
fuente]; la mezcla virgen/char es lineal en la conversión. El acero usa constantes (acero.*).
"""

import numpy as np

RHO_MIN = 50.0  # kg/m³, piso de legacy para que el lecho casi vacío no deje ρ·cp ≈ 0


def k_lecho(T, conv):
    return (1.0 - conv) * (0.12 + 1.5e-4 * (T - 298.15)) + conv * (0.08 + 2.0e-4 * (T - 298.15))


def cp_lecho(T, conv):
    return (1.0 - conv) * (1250.0 + 1.8 * (T - 273.15)) + conv * (720.0 + 2.2 * (T - 273.15))


def rho_lecho(masas):
    """masas: arreglo (n_componentes, ...) en kg/m³ (sólidos + cenizas + agua). Suma con piso de 50 kg/m³."""
    return np.maximum(np.sum(masas, axis=0), RHO_MIN)
