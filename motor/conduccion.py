"""
conduccion.py — Un paso de conducción implícita (Euler implícito) en lecho + acero.

Recibe: Malla, T (K), k (W/mK) y ρ·cp (J/m³K) por nodo, h_ext (W/m²K) por nodo axial en la cara externa,
T_horno (K), R_contacto (m²K/W, resistencia extra lecho–acero por nodo axial), dt (s).
Devuelve: (Tⁿ⁺¹, Q) con Q el calor que entró por la cara externa en el paso (J).

Ecuación (volúmenes finitos, coordenadas cilíndricas, Incropera cap. 4–5):
  ρcp·V·(Tⁿ⁺¹ − Tⁿ)/dt = Σ_vecinos G·(T_vecino − T)ⁿ⁺¹ + G_ext·(T_horno − T)ⁿ⁺¹
  G radial  = 2π r_cara·dz_k · k_cara / Δr          (k_cara = media armónica)
  G interfaz = 2π R_lecho·dz_k / (Δr_a/k_a + R_contacto)
             (el último nodo del lecho está sobre la cara r = R_lecho: del lado del lecho no hay distancia;
              el primer nodo de acero está a Δr_a de la cara. Así la resistencia no cambia al refinar la malla.
              Legacy usaba Δr_l/(2k_l) + Δr_a/(2k_a), que depende de la malla; ver lecho.R_contacto_extra.)
  G axial   = π (r_e² − r_i²) · k_cara / dz           (caras superior e inferior adiabáticas)
  G_ext     = 2π R_ext·dz_k · h_ext
En r = 0 el volumen central π(Δr/2)²·dz no tiene cara interna (equivale a L'Hôpital).
Sistema disperso de banda (ancho Nz por el orden n = i·Nz + k) resuelto con scipy.linalg.solve_banded,
que da la misma solución que spsolve y es mucho más rápido para esta malla. La energía se conserva por construcción.
"""

import numpy as np
from scipy.linalg import solve_banded


def _media_armonica(a, b):
    return 2.0 * a * b / (a + b)


def conductancias(malla, k, R_contacto):
    """Devuelve (G_r, G_z): G_r (Nr−1, Nz) entre nodos i e i+1; G_z (Nr, Nz−1) entre k y k+1. W/K."""
    m = malla
    Nl = m.n_lecho
    dr = np.diff(m.r)[:, None]
    G_r = m.A_r[1:-1, None] * m.dz_k[None, :] * _media_armonica(k[:-1], k[1:]) / dr
    if m.r.size > Nl:  # interfaz lecho–acero
        R_serie = m.dr_acero / k[Nl] + R_contacto
        G_r[Nl - 1] = m.A_r[Nl] * m.dz_k / R_serie
    area_z = np.pi * (m.r_caras[1:] ** 2 - m.r_caras[:-1] ** 2)
    G_z = area_z[:, None] * _media_armonica(k[:, :-1], k[:, 1:]) / m.dz
    return G_r, G_z


def paso_conduccion(malla, T, k, rho_cp, h_ext, T_horno, R_contacto, dt, adiabatico_z=True):
    if not adiabatico_z:
        raise NotImplementedError("Solo caras superior e inferior adiabáticas por ahora (spec §5.3)")
    m = malla
    Nr, Nz = T.shape
    G_r, G_z = conductancias(m, k, R_contacto)
    G_ext = m.A_r[-1] * m.dz_k * h_ext                     # (Nz,)
    C = rho_cp * m.V / dt                                  # (Nr, Nz)

    diag = C.copy()
    diag[:-1] += G_r
    diag[1:] += G_r
    diag[:, :-1] += G_z
    diag[:, 1:] += G_z
    diag[-1] += G_ext

    # Matriz de banda: fila u+i−j, columna j guarda A[i, j] (formato de solve_banded), u = l = Nz
    N = Nr * Nz
    ab = np.zeros((2 * Nz + 1, N))
    ab[Nz] = diag.ravel()
    acople_z = np.zeros((Nr, Nz))
    acople_z[:, :-1] = -G_z                                # A[n, n+1] = A[n+1, n] para k → k+1
    az = acople_z.ravel()[:-1]
    ab[Nz - 1, 1:] = az                                    # superdiagonal 1
    ab[Nz + 1, :-1] = az                                   # subdiagonal 1
    ar = -G_r.ravel()                                      # A[n, n+Nz] = A[n+Nz, n] para i → i+1
    ab[0, Nz:] = ar
    ab[2 * Nz, :-Nz] = ar

    b = C * T
    b[-1] += G_ext * T_horno
    T_nueva = solve_banded((Nz, Nz), ab, b.ravel(), overwrite_ab=True, check_finite=False).reshape(Nr, Nz)
    Q = float(dt * np.sum(G_ext * (T_horno - T_nueva[-1])))
    return T_nueva, Q
