"""Pruebas de calibracion/ajuste_termico.py que no necesitan correr el ajuste completo (spec de B §4 y §6)."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "calibracion"))

import ajuste_termico as aj  # noqa: E402

T_COOL = 4000       # s, último instante de la meseta


def _ensayo(sobrepaso=0.0, con_ssr=False):
    """Rampa a 500 °C, meseta y enfriamiento exponencial a 1 s, con un sobrepaso que decae al inicio de la meseta."""
    t = np.arange(0.0, 9001.0)
    t_rampa = 2850.0
    T = np.where(t < t_rampa, 25.0 + 10.0 / 60.0 * t, 500.0)
    meseta = (t >= t_rampa) & (t <= T_COOL)
    T[meseta] += sobrepaso * np.exp(-(t[meseta] - t_rampa) / 120.0)
    enf = t > T_COOL
    T[enf] = 25.0 + 475.0 * np.exp(-(t[enf] - T_COOL) / 2500.0)
    ens = {"t_s": t, "T_horno_C": T, "T_pared_C": T - 2.0, "T_centro_C": T - 10.0}
    if con_ssr:
        ens["u_ssr_pct"] = np.where(t <= T_COOL, 40.0, 0.0)
    return ens


@pytest.mark.parametrize("sobrepaso", [0.0, 8.0, 15.0])
def test_el_sobrepaso_no_adelanta_el_enfriamiento(sobrepaso):
    t = _ensayo(sobrepaso)["t_s"]
    inicio = t[aj.tramos(_ensayo(sobrepaso))["enfriamiento"]][0]
    assert T_COOL <= inicio <= T_COOL + 300.0


def test_con_ssr_el_enfriamiento_empieza_al_apagar():
    ens = _ensayo(15.0, con_ssr=True)
    assert ens["t_s"][aj.tramos(ens)["enfriamiento"]][0] == T_COOL + 1.0


def test_criterio_incluye_T_horno_de_la_corrida_1():
    c2 = {"rms": {"T_pared": 1.0, "T_centro_fuera_pirolisis": 1.0, "T_centro_pirolisis": 1.0}}
    assert aj.criterio_cierre({"rms": {"T_pared": 1.0, "T_horno": 4.9}}, c2)["T_horno_corrida1"]
    malo = aj.criterio_cierre({"rms": {"T_pared": 1.0, "T_horno": 5.1}}, c2)
    assert not malo["T_horno_corrida1"] and malo["T_pared_corrida1"]


def test_rechaza_lecho_que_no_arranca_a_T_horno():
    ens = _ensayo()
    aj.verificar_inicio({**ens, "T_centro_C": ens["T_horno_C"] + 4.0})
    with pytest.raises(ValueError, match="T_centro_C"):
        aj.verificar_inicio({**ens, "T_centro_C": ens["T_horno_C"] - 80.0})


def test_sigma_con_fijos():
    J = np.array([[1.0, 0.0], [0.0, 2.0], [1.0, 1.0]])
    cov = np.diag([0.01, 0.04])
    base = aj.sigma_con_fijos(J, cov, np.zeros((3, 1)), np.array([[1.0]]))
    assert base == pytest.approx([0.1, 0.2])
    # r = J·x + b·θ con b = J·[2, 0]: dx/dθ = [−2, 0], así que σ_θ = 0.5 suma 1 a x₀ y nada a x₁
    b = J @ np.array([2.0, 0.0])
    total = aj.sigma_con_fijos(J, cov, b[:, None], np.array([[0.25]]))
    assert total == pytest.approx([np.sqrt(0.01 + 1.0), 0.2])
