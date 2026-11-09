"""Pruebas de motor.malla: los volúmenes de control cubren exacto el lecho y el acero."""
import numpy as np
import pytest

from conftest import BASE
from motor.malla import crear_malla
from motor.parametros import cargar_parametros


@pytest.mark.parametrize("cambios", [{}, {"malla.Nr_acero": 1, "malla.Nz": 2}, {"malla.Nr_lecho": 40, "malla.Nz": 41}])
def test_volumen_total(cambios):
    p = cargar_parametros(BASE, cambios)
    m = crear_malla(p)
    R, Re, L = p["geometria.R_lecho"], p["geometria.R_lecho"] + p["geometria.e_acero"], p["geometria.L"]
    assert m.V.shape == (p["malla.Nr_lecho"] + p["malla.Nr_acero"], p["malla.Nz"])
    assert np.isclose(m.V[~m.es_acero].sum(), np.pi * R**2 * L, rtol=1e-12)
    assert np.isclose(m.V[m.es_acero].sum(), np.pi * (Re**2 - R**2) * L, rtol=1e-12)
    assert m.r[0] == 0.0 and np.isclose(m.r[m.n_lecho - 1], R) and np.isclose(m.r[-1], Re)
    assert np.isclose(m.A_r[m.n_lecho], 2 * np.pi * R) and np.isclose(m.A_r[-1], 2 * np.pi * Re) and m.A_r[0] == 0.0


def test_malla_gruesa():
    m = crear_malla(cargar_parametros(BASE, {"malla.Nr_lecho": 3}))
    assert m.n_lecho == 3
    assert np.all(m.V > 0)
    assert np.all(np.diff(m.r) > 0)
