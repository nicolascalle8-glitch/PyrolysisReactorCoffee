"""Pruebas de motor.lote y motor.graficas."""
import os

import pytest

from conftest import BASE
import motor.lote as lote_mod
from motor import Resultado, cargar_parametros, simular, simular_lote
from motor.graficas import graficar

CORTO = {"programa.rampa": 50.0, "programa.meseta": 0.0, "programa.T_final": 250.0}


def test_lote_con_falla(monkeypatch):
    original = lote_mod.simular

    def falla_si_rampa_20(p):
        if p["programa.rampa"] == 20.0:
            raise RuntimeError("falla de prueba")
        return original(p)

    monkeypatch.setattr(lote_mod, "simular", falla_si_rampa_20)
    casos = [cargar_parametros(BASE, {**CORTO, "programa.rampa": v}) for v in (50.0, 20.0, 40.0)]
    res = simular_lote(casos, procesos=1)
    assert len(res) == 3
    assert isinstance(res[0], Resultado) and isinstance(res[2], Resultado)
    assert isinstance(res[1], RuntimeError) and "falla de prueba" in str(res[1])


def test_lote_en_procesos_mantiene_orden():
    casos = [cargar_parametros(BASE, {**CORTO, "programa.T_final": T}) for T in (150.0, 250.0)]
    antes = {v: os.environ.get(v) for v in lote_mod.VARIABLES_HILOS}
    res = simular_lote(casos + [None], procesos=2)   # None: caso inválido, no tumba el lote
    assert {v: os.environ.get(v) for v in lote_mod.VARIABLES_HILOS} == antes  # el entorno queda como estaba
    assert isinstance(res[2], Exception)
    assert res[0].historia["T_set"][-1] == pytest.approx(150.0)
    assert res[1].historia["T_set"][-1] == pytest.approx(250.0)


def test_graficas(tmp_path):
    r = simular(cargar_parametros(BASE, CORTO))
    rutas = graficar(r, str(tmp_path))
    assert len(rutas) == 3
    for ruta in rutas:
        assert ruta.endswith(".png") and os.path.getsize(ruta) > 1000
