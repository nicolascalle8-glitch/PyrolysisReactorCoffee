"""Pruebas que salieron de la revisión final de la rama (2026-11-08)."""
import numpy as np
import pytest

from conftest import BASE
from motor import cargar_parametros, simular


def _rms(a, b):
    return float(np.sqrt(np.mean((a - b) ** 2)))


# Hallazgo 1 + 5: el campo de T no debe depender de la malla (la resistencia lecho–acero no puede
# cambiar al refinar) y la prueba debe poder fallar aunque los rendimientos no cambien.
@pytest.mark.lento
def test_temperatura_independiente_de_la_malla():
    p = cargar_parametros(BASE, {"programa.T_final": 300.0})  # conversión parcial
    doble = {"malla.Nr_lecho": 27, "malla.Nr_acero": 8, "malla.Nz": 39}
    a, b = simular(p), simular(p.con_cambios(doble))
    for c in ("T_centro", "T_lecho_superficie"):
        rms = _rms(a.historia[c], np.interp(a.historia["t"], b.historia["t"], b.historia[c]))
        print(f"\n{c}: RMS malla base vs ×2 = {rms:.2f} K")
        assert rms <= 2.0, f"{c}: RMS {rms:.2f} K entre malla base y ×2"
    assert abs(a.rendimientos["biochar"] - b.rendimientos["biochar"]) < 0.5


# Hallazgo 2: en modo físico el PID muestrea cada control.t_muestreo, no cada dt.
FISICO = {"horno.modo": "fisico", "horno.P_max": 5000.0, "horno.C_h": 5000.0, "control.Kd": 5.0,
          "programa.T_final": 400.0, "programa.meseta": 10.0}


@pytest.mark.lento
def test_pid_no_depende_del_dt():
    a = simular(cargar_parametros(BASE, FISICO))
    b = simular(cargar_parametros(BASE, {**FISICO, "numerica.dt_max": 1.0}))
    ma, mb = a.historia["T_termopar_control"].max(), b.historia["T_termopar_control"].max()
    assert abs(ma - mb) <= 0.3, f"sobrepico del termopar: dt_max 20 s → {ma:.2f} °C, 1 s → {mb:.2f} °C"
    seguimiento = _rms(b.historia["T_termopar_control"], b.historia["T_set"])
    assert seguimiento < 5.0, f"el PID no sigue el programa alcanzable: RMS {seguimiento:.1f} K"


# Hallazgo 3: barridos con numpy (np.arange, pandas) no deben fallar.
def test_acepta_tipos_numpy():
    p = cargar_parametros(BASE, {"malla.Nz": np.int64(20), "programa.rampa": np.float32(5.0)})
    assert p["malla.Nz"] == 20 and type(p["malla.Nz"]) is int
    assert type(p["programa.rampa"]) is float
    from motor import ErrorParametro
    with pytest.raises(ErrorParametro, match="malla.Nz"):
        cargar_parametros(BASE, {"malla.Nz": np.bool_(True)})


# Hallazgo 6: las constantes de propiedades y secado son parámetros con unidad, etiqueta y fuente.
def test_constantes_fisicas_son_parametros():
    p = cargar_parametros(BASE)
    for clave in ("lecho.k_virgen_0", "lecho.k_virgen_pend", "lecho.k_char_0", "lecho.k_char_pend",
                  "lecho.cp_virgen_0", "lecho.cp_virgen_pend", "lecho.cp_char_0", "lecho.cp_char_pend",
                  "secado.A", "secado.Ea", "secado.dH", "secado.T_umbral"):
        assert p.info(clave).fuente
        assert clave in p.supuestos() or p.info(clave).etiqueta in ("Dato", "Norma", "Calculado")
    corto = {"programa.rampa": 50.0, "programa.meseta": 0.0, "programa.T_final": 300.0}
    a = simular(cargar_parametros(BASE, corto))
    b = simular(cargar_parametros(BASE, {**corto, "lecho.k_virgen_0": 0.24}))
    assert b.historia["T_centro"][-1] > a.historia["T_centro"][-1] + 1.0
