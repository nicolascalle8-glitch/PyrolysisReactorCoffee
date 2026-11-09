"""
propiedades.py — k, cp y ρ efectivas del lecho de café.

Recibe: T en K y la conversión (0 = café virgen, 1 = char) por nodo.
Devuelve: arreglos de la misma forma.
Fuente: correlaciones lineales de legacy/simulador_interactivo.py (k_bed, cp_bed, rho_bed), [Supuesto, sin
fuente]; la mezcla virgen/char es lineal en la conversión. Los coeficientes vienen de parametros/base.toml
(lecho.k_* y lecho.cp_*); sin coeficientes se usan los de legacy. El acero usa constantes (acero.*).
  k  = (1 − X)·(k_v0 + k_vp·(T − 298.15)) + X·(k_c0 + k_cp·(T − 298.15))
  cp = (1 − X)·(cp_v0 + cp_vp·(T − 273.15)) + X·(cp_c0 + cp_cp·(T − 273.15))
"""

import numpy as np

RHO_MIN = 50.0  # kg/m³, piso de legacy para que el lecho casi vacío no deje ρ·cp ≈ 0


COEF_LEGACY = {"k_virgen_0": 0.12, "k_virgen_pend": 1.5e-4, "k_char_0": 0.08, "k_char_pend": 2.0e-4,
               "cp_virgen_0": 1250.0, "cp_virgen_pend": 1.8, "cp_char_0": 720.0, "cp_char_pend": 2.2}


def coeficientes(p):
    """Coeficientes de las correlaciones desde Parametros (claves lecho.k_* y lecho.cp_*)."""
    return {c: p[f"lecho.{c}"] for c in COEF_LEGACY}


def k_lecho(T, conv, c=COEF_LEGACY):
    return ((1.0 - conv) * (c["k_virgen_0"] + c["k_virgen_pend"] * (T - 298.15))
            + conv * (c["k_char_0"] + c["k_char_pend"] * (T - 298.15)))


def cp_lecho(T, conv, c=COEF_LEGACY):
    return ((1.0 - conv) * (c["cp_virgen_0"] + c["cp_virgen_pend"] * (T - 273.15))
            + conv * (c["cp_char_0"] + c["cp_char_pend"] * (T - 273.15)))


def rho_lecho(masas):
    """masas: arreglo (n_componentes, ...) en kg/m³ (sólidos + cenizas + agua). Suma con piso de 50 kg/m³."""
    return np.maximum(np.sum(masas, axis=0), RHO_MIN)
