"""
test_balance_masa.py
Verifica la ley de conservación de la masa en los tres modelos del proyecto:
la suma de productos + material sin reaccionar debe ser 100 % de la masa cargada,
para cualquier geometría, carga o programa térmico, y la simulación nunca debe dar NaN.

Ejecutar desde la raíz del proyecto:
    python -m pytest tests -v
"""

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "python"))

import matplotlib
matplotlib.use("Agg")

import simulador_interactivo as interactivo
from config import GeometryConfig, CoffeeBiomassConfig, HeatingRampConfig
from simulation import PyrolysisSimulator
from simulation_lumped import LumpedPyrolysisSimulator

TOLERANCIA_PCT = 0.01  # Error máximo admitido en el cierre del balance (% de la masa inicial)

# Programa corto para que las pruebas no tarden: el balance debe cerrar en cualquier instante
BASE_INTERACTIVO = dict(ramp_rate=50.0, T_final=500.0, soak_time=5.0, tau_stone=200.0,
                        R_int=0.025, e_inox=0.0025, L=0.20, w_H2O=0.08)

CASOS_INTERACTIVO = {
    "base": {},
    "canastilla_grande_R60": dict(R_int=0.060),
    "canastilla_pequena_R10": dict(R_int=0.010),
    "canastilla_alta_L40": dict(L=0.40),
    "pared_delgada_1mm": dict(e_inox=0.001),
    "pared_gruesa_6mm": dict(e_inox=0.006),
    "alta_temperatura_700": dict(T_final=700.0),
    "baja_temperatura_300": dict(T_final=300.0),
    "rampa_lenta_5": dict(ramp_rate=5.0, T_final=350.0, soak_time=0.0),
    "muy_humedo_40": dict(w_H2O=0.40),
}


@pytest.mark.parametrize("nombre", CASOS_INTERACTIVO)
def test_interactivo_cierra_balance(nombre):
    cfg = {**BASE_INTERACTIVO, **CASOS_INTERACTIVO[nombre]}
    res = interactivo.ejecutar_simulacion(cfg, verbose=False)

    assert all(np.all(np.isfinite(snap)) for snap in res["snapshots_T"]), "La temperatura divergió (NaN/inf)"

    total = (res["pct_char"] + res["pct_vol"] + res["pct_gnc"] + res["pct_water"]
             + res["pct_biomass"] + res["pct_water_remaining"])
    assert np.all(np.abs(total - 100.0) < TOLERANCIA_PCT), f"Balance no cierra: {total.min():.4f} a {total.max():.4f} %"
    for clave in ("pct_char", "pct_vol", "pct_gnc", "pct_water", "pct_biomass", "pct_water_remaining"):
        assert np.all(res[clave] >= 0.0), f"{clave} negativo"


def _sim_python(Sim, geo=None, bio=None, ramp=None):
    sim = Sim(geo_cfg=GeometryConfig(**(geo or {})),
              bio_cfg=CoffeeBiomassConfig(**(bio or {})),
              ramp_cfg=HeatingRampConfig(**(ramp or {})))
    return sim.run(verbose=False)


CASOS_PYTHON = {
    "canastilla_grande_R60": dict(geo=dict(R_inner=0.06, R_outer_tube=0.062),
                                  ramp=dict(ramp_rate_C_per_min=50.0, soak_time_min=5.0)),
    "alta_temperatura_700": dict(ramp=dict(ramp_rate_C_per_min=50.0, T_final_C=700.0, soak_time_min=5.0)),
}


@pytest.mark.parametrize("Sim", [PyrolysisSimulator, LumpedPyrolysisSimulator], ids=["multicomponente", "lumped"])
@pytest.mark.parametrize("nombre", CASOS_PYTHON)
def test_modelos_python_cierran_balance(Sim, nombre):
    res = _sim_python(Sim, **CASOS_PYTHON[nombre])
    fy = res["final_yields"]
    assert all(np.isfinite(v) and v >= 0.0 for v in fy.values()), fy
    assert abs(sum(fy.values()) - 100.0) < TOLERANCIA_PCT, f"Balance = {sum(fy.values()):.4f} %: {fy}"
