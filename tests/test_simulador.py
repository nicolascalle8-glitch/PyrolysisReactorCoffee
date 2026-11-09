"""
Pruebas de motor.simulador: caso base contra legacy (prueba 2), balances en 16 casos (prueba 3),
temperatura baja, divergencia y exportación.
"""
import json

import numpy as np
import pytest

from conftest import BASE
import motor.simulador as sim_mod
from motor import ErrorSimulacion, cargar_parametros, simular

LEGACY_REF = {"biochar": 24.93, "bio_oil": 38.94, "gas": 20.23, "agua": 15.90}  # spec §9, prueba 2

# Casos de legacy/tests/test_balance_masa.py (CASOS_INTERACTIVO) traducidos a cambios del TOML,
# sobre el mismo programa corto de esas pruebas (50 °C/min, meseta 5 min).
CORTO = {"programa.rampa": 50.0, "programa.meseta": 5.0}
CASOS = {
    "base": {},
    "canastilla_grande_R60": {"geometria.R_lecho": 0.060},
    "canastilla_pequena_R10": {"geometria.R_lecho": 0.010},
    "canastilla_alta_L40": {"geometria.L": 0.40},
    "pared_delgada_1mm": {"geometria.e_acero": 0.001},
    "pared_gruesa_6mm": {"geometria.e_acero": 0.006},
    "alta_temperatura_700": {"programa.T_final": 700.0},
    "baja_temperatura_300": {"programa.T_final": 300.0},
    "rampa_lenta_5": {"programa.rampa": 5.0, "programa.T_final": 350.0, "programa.meseta": 0.0},
    "muy_humedo_40": {"lecho.humedad": 40.0},
    "rac_regimen_rapida": {"cinetica.regimen": "rapida"},
    "acetilos_hce_04": {"cinetica.acetilos_hce": 0.4},   # reemplaza 'cinetica_original' (no existe en el motor)
    # 4 de legacy/python/simulation.py (programa completo)
    "sim_rampa_5": {"programa.rampa": 5.0, "programa.meseta": 30.0},
    "sim_rampa_20": {"programa.rampa": 20.0, "programa.meseta": 30.0},
    "sim_R60": {"geometria.R_lecho": 0.060, "programa.rampa": 10.0, "programa.meseta": 30.0},
    "sim_humedad_40": {"lecho.humedad": 40.0, "programa.rampa": 10.0, "programa.meseta": 30.0},
}


@pytest.fixture(scope="module")
def base():
    # Tarea 8 (dt fijo): con dt = 2 s el paso dividido corre ~0.3 min la exoterma del centro (RMS 13 K);
    # con dt → 0 converge a RMS ≈ 6 K. Se usa dt0 = 0.25 s hasta que la Tarea 9 meta el dt adaptativo.
    return simular(cargar_parametros(BASE, {"numerica.dt0": 0.25}))


@pytest.fixture(scope="module")
def legacy_base():
    import simulador_interactivo as legacy
    cfg = dict(cinetica="rac", regimen="lecho_fijo", ramp_rate=10.0, T_final=500.0, soak_time=30.0,
               tau_stone=200.0, R_int=0.025, e_inox=0.0025, L=0.20, w_H2O=0.08)
    return legacy.ejecutar_simulacion(cfg, verbose=False)


@pytest.mark.lento
def test_caso_base_vs_legacy(base, legacy_base):
    r = base.rendimientos
    for clave, ref in LEGACY_REF.items():
        assert abs(r[clave] - ref) <= 1.0, f"{clave}: motor {r[clave]:.2f} vs legacy {ref:.2f}"
    # la referencia del spec debe ser lo que legacy da hoy
    assert legacy_base["pct_char"][-1] == pytest.approx(LEGACY_REF["biochar"], abs=0.01)
    h = base.historia
    T_leg = np.interp(h["t"], legacy_base["time"], legacy_base["T_center"])
    rms = float(np.sqrt(np.mean((h["T_centro"] - T_leg) ** 2)))
    assert rms <= 10.0, f"RMS T_centro = {rms:.2f} K"


@pytest.mark.parametrize("nombre", CASOS)
def test_balances_16_casos(nombre):
    cambios = dict(CASOS[nombre]) if nombre.startswith("sim_") else {**CORTO, **CASOS[nombre]}
    r = simular(cargar_parametros(BASE, cambios))
    assert r.balance_masa <= 0.01, f"masa {r.balance_masa:.2e} %"
    assert r.balance_energia <= 0.5, f"energía {r.balance_energia:.2e} %"
    assert r.confiable
    assert all(np.isfinite(v) and v >= 0 for v in r.rendimientos.values())


def test_temperatura_baja():
    p = cargar_parametros(BASE, {"programa.T_final": 150.0})
    r = simular(p)
    # La carga húmeda trae agua y cenizas: la biomasa original es (1 − w)(1 − cenizas) ≈ 89.9 % de la carga.
    organico = 100.0 * (1 - p["lecho.humedad"] / 100) * (1 - p["lecho.cenizas"] / 100)
    assert r.rendimientos["sin_reaccionar"] > 0.9 * organico
    assert r.confiable


def test_humedad_cero():
    r = simular(cargar_parametros(BASE, {**CORTO, "lecho.humedad": 0.0}))
    assert r.rendimientos["humedad_restante"] == 0.0 and r.confiable


def test_modo_fisico_onoff_corre():
    r = simular(cargar_parametros(BASE, {**CORTO, "horno.modo": "fisico", "control.modo_ssr": "onoff"}))
    assert r.confiable
    assert np.all(np.isfinite(r.historia["T_horno"]))


def test_divergencia(monkeypatch):
    original = sim_mod.paso_conduccion

    def con_nan(*args, **kw):
        T, Q = original(*args, **kw)
        T[3, 2] = np.nan
        return T, Q

    monkeypatch.setattr(sim_mod, "paso_conduccion", con_nan)
    with pytest.raises(ErrorSimulacion) as e:
        simular(cargar_parametros(BASE, CORTO))
    assert "t =" in str(e.value) and "nodo" in str(e.value)


def test_a_csv(tmp_path):
    r = simular(cargar_parametros(BASE, {**CORTO, "programa.T_final": 200.0}))
    r.a_csv(str(tmp_path / "caso"))
    for nombre in ("historia.csv", "rendimientos.json", "mapas.npz"):
        assert (tmp_path / "caso" / nombre).exists()
    datos = json.loads((tmp_path / "caso" / "rendimientos.json").read_text(encoding="utf-8"))
    assert set(datos["rendimientos"]) == {"biochar", "bio_oil", "gas", "agua", "sin_reaccionar", "humedad_restante"}
    assert set(r.historia) >= {"t", "T_horno", "T_set", "T_pared", "T_centro", "T_termopar_control", "masa_residual"}
