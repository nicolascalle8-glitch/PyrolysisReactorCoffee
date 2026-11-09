"""Pruebas de horno.modo = "medido": T_horno(t) leído de un CSV (decisión 11, spec de B §5)."""
import csv

import numpy as np
import pytest

from conftest import BASE
from motor import cargar_parametros, simular
from motor.horno import HornoMedido
from motor.parametros import ErrorParametro

RAPIDO = {"programa.T_final": 160.0, "programa.meseta": 5.0, "malla.Nr_lecho": 6, "malla.Nr_acero": 2, "malla.Nz": 3}


def _csv(ruta, t, T_h, extra=None):
    extra = extra or {}
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["t_s", "T_horno_C", *extra])
        for i, a in enumerate(t):
            w.writerow([a, T_h[i], *[v[i] for v in extra.values()]])
    return str(ruta)


def test_lee_y_interpola(tmp_path):
    h = HornoMedido(_csv(tmp_path / "a.csv", [10.0, 20.0, 30.0], [100.0, 200.0, 150.0],
                         {"T_set_C": [90.0, 95.0, 100.0], "u_ssr_pct": [0.0, 50.0, 100.0]}))
    assert h.t_fin == 20.0                                  # el tiempo cuenta desde el primer registro
    assert h.temperatura(5.0) == pytest.approx(150.0 + 273.15)
    assert h.setpoint(10.0) == pytest.approx(95.0 + 273.15)
    assert h.salida(15.0) == pytest.approx(0.75)


def test_columnas_opcionales_quedan_nan(tmp_path):
    h = HornoMedido(_csv(tmp_path / "a.csv", [0.0, 1.0], [25.0, 26.0]))
    assert np.isnan(h.setpoint(0.5)) and np.isnan(h.salida(0.5))


@pytest.mark.parametrize("t, T, patron", [
    ([0.0, 0.0, 1.0], [25.0, 26.0, 27.0], "creciente"),
    ([0.0, 1.0], ["", 26.0], "fila 2"),
    ([0.0], [25.0], "2 registros"),
])
def test_csv_invalido(tmp_path, t, T, patron):
    with pytest.raises(ErrorParametro, match=patron):
        HornoMedido(_csv(tmp_path / "a.csv", t, T))


def test_columna_faltante(tmp_path):
    ruta = tmp_path / "a.csv"
    ruta.write_text("t_s,T_pared_C\n0,25\n1,26\n", encoding="utf-8")
    with pytest.raises(ErrorParametro, match="T_horno_C"):
        HornoMedido(str(ruta))


def test_archivo_no_existe(tmp_path):
    with pytest.raises(ErrorParametro, match="horno.archivo_medido"):
        HornoMedido(str(tmp_path / "no_existe.csv"))


def test_modo_medido_exige_archivo():
    with pytest.raises(ErrorParametro, match="horno.archivo_medido"):
        cargar_parametros(BASE, {"horno.modo": "medido"})


def test_el_archivo_debe_ser_texto():
    with pytest.raises(ErrorParametro, match="horno.archivo_medido"):
        cargar_parametros(BASE, {"horno.archivo_medido": 3.0})


def test_sigue_el_csv_y_termina_con_el(tmp_path):
    t = np.arange(0.0, 1201.0)
    T_h = 25.0 + 0.15 * t
    ruta = _csv(tmp_path / "a.csv", t, T_h, {"T_set_C": T_h + 5.0, "u_ssr_pct": np.full(t.size, 40.0)})
    h = simular(cargar_parametros(BASE, {**RAPIDO, "horno.modo": "medido", "horno.archivo_medido": ruta})).historia
    assert h["t"][-1] * 60.0 == pytest.approx(1200.0)
    assert h["T_horno"] == pytest.approx(np.interp(h["t"] * 60.0, t, T_h), abs=1e-6)
    assert h["T_set"] == pytest.approx(h["T_horno"] + 5.0, abs=1e-6)
    assert h["u"] == pytest.approx(0.4)
    assert h["T_pared"][0] == pytest.approx(25.0)           # arranca a la T_horno medida
    assert np.all(h["T_centro"] <= h["T_pared"] + 1e-9)


def test_reproduce_el_modo_ideal(tmp_path):
    """Con T_horno(t) del modo ideal como entrada, el tubo y el lecho deben dar lo mismo que en modo ideal."""
    base = {**RAPIDO, "numerica.t_registro": 1.0, "numerica.dt_max": 2.0}
    ideal = simular(cargar_parametros(BASE, base)).historia
    ruta = _csv(tmp_path / "a.csv", ideal["t"] * 60.0, ideal["T_horno"])
    medido = simular(cargar_parametros(BASE, {**base, "horno.modo": "medido", "horno.archivo_medido": ruta})).historia
    ts = medido["t"]
    for c in ("T_pared", "T_centro"):
        assert np.max(np.abs(medido[c] - np.interp(ts, ideal["t"], ideal[c]))) < 0.05, c
