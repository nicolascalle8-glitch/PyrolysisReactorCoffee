"""Pruebas de motor.cinetica: RAC + TGL (traídas de legacy/tests/test_cinetica_rac.py) y secado."""
import numpy as np
import pytest

from motor.cinetica import REGIMENES, CineticaRAC, paso_secado


@pytest.mark.parametrize("regimen", REGIMENES)
def test_reacciones_balanceadas_atomo_a_atomo(regimen):
    kin = CineticaRAC(regimen)
    desbalanceadas = {n: e for n, e in kin.error_atomico.items() if e > 1e-6}
    assert not desbalanceadas, f"Reacciones con desbalance C/H/O: {desbalanceadas}"
    assert np.allclose(kin.Y.sum(axis=0), 1.0)


@pytest.mark.parametrize("dt", [1e-3, 1.0, 1e3])
@pytest.mark.parametrize("T", [500.0, 700.0, 1100.0])
def test_masa_se_conserva_para_cualquier_paso(T, dt):
    kin = CineticaRAC("lecho_fijo")
    m = kin.estado_inicial(1.0, forma=(3, 4))
    total = {"gas": 0.0, "agua_rxn": 0.0, "bio_oil": 0.0}
    for _ in range(20):
        _, lib = kin.paso(m, np.full((3, 4), T), dt)
        for f in total:
            total[f] += lib[f]
    assert np.all(m >= 0.0)
    assert np.allclose(kin.masa_solida(m) + sum(total.values()), 1.0, atol=1e-12)


def test_igual_a_legacy():
    """La química no cambió: mismo resultado que legacy/python/kinetics_rac.py."""
    import kinetics_rac as viejo  # legacy/conftest.py no aplica aquí: se agrega la ruta a mano
    for regimen in REGIMENES:
        a, b = CineticaRAC(regimen), viejo.CineticaRAC(regimen)
        ma, mb = a.estado_inicial(300.0, forma=(2,)), b.estado_inicial(300.0, forma=(2,))
        qa, la = a.paso(ma, np.array([600.0, 750.0]), 10.0)
        qb, lb = b.paso(mb, np.array([600.0, 750.0]), 10.0)
        assert np.array_equal(ma, mb) and np.array_equal(qa, qb)
        assert all(np.array_equal(la[f], lb[f]) for f in la)


def test_secado_humedad_cero():
    agua, evap, q = paso_secado(np.zeros(5), np.full(5, 450.0), 2.0)
    assert np.all(np.isfinite(agua)) and np.all(np.isfinite(q))
    assert np.all(evap == 0.0) and np.all(q == 0.0)


def test_secado_conserva():
    agua = np.linspace(0.0, 40.0, 7)
    T = np.linspace(300.0, 600.0, 7)
    nueva, evap, q = paso_secado(agua, T, 3.0)
    assert np.allclose(nueva + evap, agua, rtol=1e-12, atol=0)
    assert np.all(evap[T <= 350.0] == 0.0)            # umbral de legacy
    assert np.all(q <= 0.0)                           # evaporar enfría
    assert np.allclose(q, evap / 3.0 * -2.26e6)
