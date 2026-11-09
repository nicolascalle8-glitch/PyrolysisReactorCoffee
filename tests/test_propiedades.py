"""Pruebas de motor.propiedades: mismas correlaciones que legacy/simulador_interactivo.py (líneas 233–235)."""
import itertools

import numpy as np

from motor.propiedades import cp_lecho, k_lecho, rho_lecho


def _legacy(T, conv):
    k = (1.0 - conv) * (0.12 + 1.5e-4 * (T - 298.15)) + conv * (0.08 + 2.0e-4 * (T - 298.15))
    cp = (1.0 - conv) * (1250.0 + 1.8 * (T - 273.15)) + conv * (720.0 + 2.2 * (T - 273.15))
    return k, cp


def test_igual_a_legacy():
    pares = list(itertools.product([300.0, 500.0, 800.0], [0.0, 0.5, 1.0]))
    T = np.array([a for a, _ in pares])
    c = np.array([b for _, b in pares])
    k, cp = _legacy(T, c)
    assert np.allclose(k_lecho(T, c), k, rtol=1e-12, atol=0)
    assert np.allclose(cp_lecho(T, c), cp, rtol=1e-12, atol=0)


def test_rho_con_piso():
    masas = np.array([[300.0, 20.0, 1.0], [10.0, 5.0, 0.0]])  # filas = nodos, columnas = componentes
    assert np.allclose(rho_lecho(masas.T), [321.0, 50.0])
