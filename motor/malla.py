"""
malla.py — Nodos r–z del lecho y del tubo de acero, con sus volúmenes de control.

Recibe: Parametros (geometria.*, malla.*).
Devuelve: Malla.
Nodos radiales como legacy: lecho en linspace(0, R_lecho, Nr_lecho) (nodo 0 en r = 0, último en la
pared de la canastilla) y acero en R_lecho + j·e/Nr_acero, j = 1..Nr_acero (último en la cara externa).
Las caras de los volúmenes están a medio camino entre nodos; las caras extremas son r = 0, r = R_lecho
(interfaz lecho–acero) y r = R_ext. Axial: nodos en linspace(0, L, Nz) con medio volumen en los extremos.
Así los volúmenes suman exacto π R² L (lecho) y π (R_ext² − R²) L (acero).
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class Malla:
    r: np.ndarray            # (Nr_total,) radio de cada nodo, m
    z: np.ndarray            # (Nz,) altura de cada nodo, m
    n_lecho: int             # nodos del lecho (los primeros n_lecho de r)
    es_acero: np.ndarray     # (Nr_total, Nz) bool
    V: np.ndarray            # (Nr_total, Nz) volumen de control, m³
    A_r: np.ndarray          # (Nr_total + 1,) área de cara radial por unidad de dz = 2π r_cara, m
    r_caras: np.ndarray      # (Nr_total + 1,) radio de cada cara radial, m
    dz: float                # distancia entre nodos axiales, m
    dz_k: np.ndarray         # (Nz,) espesor axial de cada volumen (dz, o dz/2 en los extremos), m
    dr_lecho: float
    dr_acero: float
    R_lecho: float
    R_ext: float

    @property
    def forma(self):
        return self.V.shape


def crear_malla(p, con_acero: bool = True) -> Malla:
    """con_acero=False deja solo el lecho (cara externa en R_lecho); lo usa la prueba de Bessel."""
    R = p["geometria.R_lecho"]
    e = p["geometria.e_acero"] if con_acero else 0.0
    L = p["geometria.L"]
    Nl, Nz = p["malla.Nr_lecho"], p["malla.Nz"]
    Na = p["malla.Nr_acero"] if con_acero else 0
    R_ext = R + e

    r_l = np.linspace(0.0, R, Nl)
    dr_a = e / Na if Na else 0.0
    r_a = R + dr_a * np.arange(1, Na + 1)
    if Na:
        r_a[-1] = R_ext
    r = np.concatenate([r_l, r_a])

    caras = np.empty(Nl + Na + 1)
    caras[0] = 0.0
    caras[1:Nl] = 0.5 * (r_l[1:] + r_l[:-1])
    caras[Nl] = R
    if Na:
        caras[Nl + 1:Nl + Na] = 0.5 * (r_a[1:] + r_a[:-1])
    caras[-1] = R_ext

    z = np.linspace(0.0, L, Nz)
    dz = L / (Nz - 1)
    dz_k = np.full(Nz, dz)
    dz_k[0] = dz_k[-1] = dz / 2.0

    area = np.pi * (caras[1:] ** 2 - caras[:-1] ** 2)
    V = area[:, None] * dz_k[None, :]
    es_acero = np.zeros((Nl + Na, Nz), dtype=bool)
    es_acero[Nl:] = True
    return Malla(r=r, z=z, n_lecho=Nl, es_acero=es_acero, V=V, A_r=2.0 * np.pi * caras, r_caras=caras,
                 dz=dz, dz_k=dz_k, dr_lecho=r_l[1] - r_l[0], dr_acero=dr_a, R_lecho=R, R_ext=R_ext)
