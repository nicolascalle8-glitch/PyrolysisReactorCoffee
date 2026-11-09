"""
Pruebas de motor.conduccion (prueba 1 del spec y conservación de energía).

Bessel: cilindro infinito de lecho homogéneo, T inicial uniforme, pared a T fija (h_ext muy grande).
Solución exacta: θ = (T − T∞)/(T0 − T∞) = Σ 2/(λn J1(λn)) J0(λn r/R) exp(−λn² Fo), con λn raíces de J0
(Incropera, Fundamentals of Heat and Mass Transfer, sección 5.6, conducción transitoria en cilindro;
edición a confirmar con el libro de Nick).
"""
import numpy as np
import pytest
from scipy.special import j0, j1, jn_zeros

from conftest import BASE
from motor.conduccion import paso_conduccion
from motor.malla import crear_malla
from motor.parametros import cargar_parametros

K, RHO_CP = 0.1, 5.0e5          # lecho homogéneo de la prueba (W/mK, J/m³K)
T0, T_INF = 300.0, 400.0


def _theta_exacta(r, R, Fo, n=50):
    lam = jn_zeros(0, n)
    return np.sum(2.0 / (lam * j1(lam)) * j0(np.outer(r / R, lam)) * np.exp(-lam**2 * Fo), axis=1)


def _campos(m, k=K, rho_cp=RHO_CP):
    return np.full(m.forma, k), np.full(m.forma, rho_cp)


@pytest.mark.parametrize("Fo", [0.05, 0.2, 0.5])
def test_bessel(Fo):
    p = cargar_parametros(BASE, {"malla.Nr_lecho": 30, "malla.Nz": 3})
    m = crear_malla(p, con_acero=False)
    R = m.R_lecho
    k, rc = _campos(m)
    t_fin = Fo * R**2 * RHO_CP / K
    n_pasos = 400
    dt = t_fin / n_pasos
    T = np.full(m.forma, T0)
    h = np.full(m.forma[1], 1e9)
    Rc = np.zeros(m.forma[1])
    for _ in range(n_pasos):
        T, _ = paso_conduccion(m, T, k, rc, h, T_INF, Rc, dt)
    exacta = T_INF + (T0 - T_INF) * _theta_exacta(m.r, R, Fo)
    error = np.max(np.abs(T[:, 1] - exacta))
    assert error <= 0.01 * (T_INF - T0), f"Fo = {Fo}: error máximo {error:.3f} K"


def _caso_con_acero():
    p = cargar_parametros(BASE)
    m = crear_malla(p)
    k, rc = _campos(m)
    k[m.es_acero] = p["acero.k"]
    rc[m.es_acero] = p["acero.rho"] * p["acero.cp"]
    T = 300.0 + 50.0 * np.random.default_rng(1).random(m.forma)
    return p, m, k, rc, T


def test_energia_adiabatica():
    p, m, k, rc, T = _caso_con_acero()
    E0 = np.sum(rc * m.V * T)
    Rc = np.full(m.forma[1], p["lecho.R_contacto_extra"])
    for _ in range(20):
        T, Q = paso_conduccion(m, T, k, rc, np.zeros(m.forma[1]), 900.0, Rc, 5.0)
        assert Q == 0.0
    assert np.isclose(np.sum(rc * m.V * T), E0, rtol=1e-10)


def test_flujo_interfaz():
    p, m, k, rc, T = _caso_con_acero()
    E0 = np.sum(rc * m.V * T)
    Rc = np.full(m.forma[1], p["lecho.R_contacto_extra"])
    Q_total = 0.0
    for _ in range(20):
        T, Q = paso_conduccion(m, T, k, rc, np.full(m.forma[1], 90.0), 900.0, Rc, 5.0)
        Q_total += Q
    assert Q_total > 0
    assert np.isclose(np.sum(rc * m.V * T) - E0, Q_total, rtol=1e-9)
