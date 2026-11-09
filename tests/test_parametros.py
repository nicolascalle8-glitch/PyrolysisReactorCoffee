"""Pruebas de motor.parametros: carga, cambios y validación (prueba 8 del spec)."""
import pytest

from conftest import BASE
from motor.parametros import ErrorParametro, cargar_parametros


def test_carga_base():
    p = cargar_parametros(BASE)
    assert p["programa.rampa"] == 10.0
    assert p.info("lecho.cenizas").etiqueta == "Dato"
    assert "lecho.hemicelulosa" in p.supuestos() and "lecho.celulosa" not in p.supuestos()


def test_cambio_valido():
    assert cargar_parametros(BASE, {"control.Kp": 4.0})["control.Kp"] == 4.0


def test_cambio_no_toca_el_original():
    p = cargar_parametros(BASE)
    q = p.con_cambios({"programa.rampa": 20.0})
    assert q["programa.rampa"] == 20.0 and p["programa.rampa"] == 10.0


def test_clave_inexistente():
    with pytest.raises(ErrorParametro, match="control.kp"):
        cargar_parametros(BASE, {"control.kp": 4})


def test_fuera_de_rango():
    with pytest.raises(ErrorParametro, match="malla.Nr_lecho"):
        cargar_parametros(BASE, {"malla.Nr_lecho": 1})
    assert cargar_parametros(BASE, {"malla.Nr_lecho": 3})["malla.Nr_lecho"] == 3


def test_unidad_mala(tmp_path):
    texto = open(BASE, encoding="utf-8").read()
    malo = texto.replace('unidad   = "°C/min"', 'unidad   = "K/s"', 1)
    assert malo != texto
    ruta = tmp_path / "malo.toml"
    ruta.write_text(malo, encoding="utf-8")
    with pytest.raises(ErrorParametro, match="unidad") as e:
        cargar_parametros(str(ruta))
    assert "programa.rampa" in str(e.value)


def test_composicion_no_suma():
    with pytest.raises(ErrorParametro, match="composición"):
        cargar_parametros(BASE, {"lecho.celulosa": 50})


def test_opcion_invalida():
    with pytest.raises(ErrorParametro, match="horno.modo"):
        cargar_parametros(BASE, {"horno.modo": "magico"})


def test_composicion_con_fuente_ballesteros():
    """Composición aprobada por Nick (2026-11-08) tras la REVISIÓN de la investigación de Gemini."""
    p = cargar_parametros(BASE)
    assert p["lecho.celulosa"] == 12.4 and p["lecho.lignina"] == 23.9
    for c in ("celulosa", "lignina"):
        assert p.info(f"lecho.{c}").etiqueta == "Dato" and "Ballesteros" in p.info(f"lecho.{c}").fuente
    assert p.info("lecho.lipidos").etiqueta == "Por validar"
    assert p.info("lecho.rho_aparente").etiqueta == "Por validar"
