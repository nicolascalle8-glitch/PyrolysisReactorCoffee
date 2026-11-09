"""
Validación de la cinética RAC + TGL contra datos de café (traída de legacy/tests/test_cinetica_rac.py).

Muestra 0D térmicamente delgada (como un TGA): T impuesta (rampa + meseta), solo cambia la química.
  - Bok et al. (2012), Energy 47(1) 17-24: pirólisis rápida, máximo de bio-oil 54.85 % a 823 K.
  - Inayat et al. (2022), Environ. Clim. Technol. 26(1) 720-729: lecho fijo, ~0.9 °C/min hasta 300 °C,
    30 min, café seco → bio-oil ~22.5 %.
Tolerancias iguales a las de legacy (±3 y ±5 puntos).
"""
import pytest

from motor.cinetica.tga import simular_0d

BOK = dict(regimen="rapida", rampa=60000.0, T_final=550.0, meseta=1.0, humedad=0.0)
INAYAT = dict(regimen="lecho_fijo", rampa=0.9, T_final=300.0, meseta=30.0, humedad=0.0)


def test_validacion_pirolisis_rapida_bok_2012():
    r = simular_0d(**BOK)
    assert abs(r["bio_oil"] + r["agua"] - 54.85) < 3.0


def test_validacion_lecho_fijo_inayat_2022_liquido():
    r = simular_0d(**INAYAT)
    assert abs(r["bio_oil"] + r["agua"] - 22.5) < 5.0


def test_rapida_produce_mas_liquido_que_lenta():
    rapida = simular_0d(**BOK)
    lenta = simular_0d(**dict(BOK, rampa=10.0, meseta=30.0, regimen="lecho_fijo"))
    assert rapida["bio_oil"] > lenta["bio_oil"]
    assert rapida["biochar"] < lenta["biochar"]


def test_igual_a_legacy_0d():
    """Mismos números que legacy/python/comparar_cineticas.simular_rac (corrido el 2026-11-08)."""
    assert simular_0d(**INAYAT)["bio_oil"] == pytest.approx(18.250, abs=0.01)
    assert simular_0d(**BOK)["bio_oil"] == pytest.approx(50.485, abs=0.01)


def test_cierra_masa():
    r = simular_0d(**dict(INAYAT, humedad=8.0))
    assert sum(r.values()) == pytest.approx(100.0, abs=1e-9)
