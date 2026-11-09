"""Pruebas 4, 5 y 6 del spec: convergencia en dt, convergencia de malla y velocidad."""
import os
import platform
import time

import pytest

from conftest import BASE
from motor import cargar_parametros, simular

CLAVES = ("biochar", "bio_oil", "gas", "agua")


def _max_dif(a, b):
    return max(abs(a.rendimientos[c] - b.rendimientos[c]) for c in CLAVES)


@pytest.mark.lento
def test_dt_mitad():
    p = cargar_parametros(BASE)
    mitad = {c: p[c] / 2 for c in ("numerica.dt0", "numerica.dt_max", "numerica.dT_max", "numerica.dX_max")}
    a, b = simular(p), simular(p.con_cambios(mitad))
    d = _max_dif(a, b)
    print(f"\nprueba 4: dt/2 cambia los rendimientos {d:.3f} pp (pasos {a.pasos} → {b.pasos})")
    assert d < 0.2


@pytest.mark.lento
def test_malla_doble():
    p = cargar_parametros(BASE)
    doble = {"malla.Nr_lecho": 2 * p["malla.Nr_lecho"] - 1, "malla.Nr_acero": 2 * p["malla.Nr_acero"],
             "malla.Nz": 2 * p["malla.Nz"] - 1}
    a, b = simular(p), simular(p.con_cambios(doble))
    rango = {c: (round(a.rendimientos[c], 2), round(b.rendimientos[c], 2)) for c in CLAVES}
    d = _max_dif(a, b)
    print(f"\nprueba 5: malla ×2 cambia los rendimientos {d:.3f} pp; rango (base, ×2): {rango}")
    assert d < 0.5


def test_rampa_rapida_meseta_cero():
    t0 = time.perf_counter()
    r = simular(cargar_parametros(BASE, {"programa.rampa": 50.0, "programa.meseta": 0.0}))
    assert time.perf_counter() - t0 < 60.0
    assert r.confiable
    assert r.historia["t"][-1] == pytest.approx((500 - 25) / 50.0)


def test_dt_adaptativo_respeta_limites():
    p = cargar_parametros(BASE, {"programa.rampa": 50.0, "programa.meseta": 5.0})
    r = simular(p)
    assert r.dT_max_paso <= p["numerica.dT_max"] + 1e-9 or r.pasos_en_dt_min > 0


@pytest.mark.lento
def test_velocidad():
    p = cargar_parametros(BASE)
    t0 = time.perf_counter()
    simular(p)
    dur = time.perf_counter() - t0
    maquina = f"{platform.machine()} {platform.processor() or platform.system()}, {os.cpu_count()} núcleos"
    assert dur <= 5.0, f"caso base tardó {dur:.1f} s en {maquina} (meta ≤ 5 s en un núcleo del M4; legacy 25 s)"
