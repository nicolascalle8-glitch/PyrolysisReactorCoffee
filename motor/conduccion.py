"""
conduccion.py — Un paso de conducción implícita (Euler implícito) en lecho + acero.

Recibe: Malla, T (K), k (W/mK) y ρ·cp (J/m³K) por nodo, h_ext (W/m²K) por nodo axial en la cara externa,
T_horno (K), R_contacto (m²K/W, resistencia extra lecho–acero por nodo axial), dt (s).
Devuelve: (Tⁿ⁺¹, Q) con Q el calor que entró por la cara externa en el paso (J).

Ecuación (volúmenes finitos, coordenadas cilíndricas, Incropera cap. 4–5):
  ρcp·V·(Tⁿ⁺¹ − Tⁿ)/dt = Σ_vecinos G·(T_vecino − T)ⁿ⁺¹ + G_ext·(T_horno − T)ⁿ⁺¹
  G radial  = 2π r_cara·dz_k · k_cara / Δr          (k_cara = media armónica)
  G interfaz = 2π R_lecho·dz_k / (Δr_l/(2k_l) + Δr_a/(2k_a) + R_contacto)   (como legacy)
  G axial   = π (r_e² − r_i²) · k_cara / dz           (caras superior e inferior adiabáticas)
  G_ext     = 2π R_ext·dz_k · h_ext
En r = 0 el volumen central π(Δr/2)²·dz no tiene cara interna (equivale a L'Hôpital).
Sistema disperso (scipy.sparse) resuelto con spsolve. La energía se conserva por construcción.
"""

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve


def _media_armonica(a, b):
    return 2.0 * a * b / (a + b)


def conductancias(malla, k, R_contacto):
    """Devuelve (G_r, G_z): G_r (Nr−1, Nz) entre nodos i e i+1; G_z (Nr, Nz−1) entre k y k+1. W/K."""
    m = malla
    Nl = m.n_lecho
    dr = np.diff(m.r)[:, None]
    G_r = m.A_r[1:-1, None] * m.dz_k[None, :] * _media_armonica(k[:-1], k[1:]) / dr
    if m.r.size > Nl:  # interfaz lecho–acero
        R_serie = m.dr_lecho / (2.0 * k[Nl - 1]) + m.dr_acero / (2.0 * k[Nl]) + R_contacto
        G_r[Nl - 1] = m.A_r[Nl] * m.dz_k / R_serie
    area_z = np.pi * (m.r_caras[1:] ** 2 - m.r_caras[:-1] ** 2)
    G_z = area_z[:, None] * _media_armonica(k[:, :-1], k[:, 1:]) / m.dz
    return G_r, G_z


def paso_conduccion(malla, T, k, rho_cp, h_ext, T_horno, R_contacto, dt, adiabatico_z=True):
    m = malla
    Nr, Nz = T.shape
    G_r, G_z = conductancias(m, k, R_contacto)
    if not adiabatico_z:
        raise NotImplementedError("Solo caras superior e inferior adiabáticas por ahora (spec §5.3)")
    G_ext = m.A_r[-1] * m.dz_k * h_ext                     # (Nz,)
    C = rho_cp * m.V / dt                                  # (Nr, Nz)

    idx = np.arange(Nr * Nz).reshape(Nr, Nz)
    diag = C.copy()
    filas, cols, vals = [], [], []

    def acoplar(a, b, G):
        nonlocal diag
        filas.extend([a.ravel(), b.ravel()])
        cols.extend([b.ravel(), a.ravel()])
        vals.extend([-G.ravel(), -G.ravel()])

    acoplar(idx[:-1], idx[1:], G_r)
    diag[:-1] += G_r
    diag[1:] += G_r
    acoplar(idx[:, :-1], idx[:, 1:], G_z)
    diag[:, :-1] += G_z
    diag[:, 1:] += G_z
    diag[-1] += G_ext

    b = C * T
    b[-1] += G_ext * T_horno
    filas.append(idx.ravel()); cols.append(idx.ravel()); vals.append(diag.ravel())
    A = csc_matrix((np.concatenate(vals), (np.concatenate(filas), np.concatenate(cols))), shape=(Nr * Nz,) * 2)
    T_nueva = spsolve(A, b.ravel()).reshape(Nr, Nz)
    Q = float(dt * np.sum(G_ext * (T_horno - T_nueva[-1])))
    return T_nueva, Q
