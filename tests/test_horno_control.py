"""Pruebas de motor.horno, motor.control y motor.termopar."""
import numpy as np
import pytest

from conftest import BASE
from motor.control import PID, SSR, Programa
from motor.horno import HornoFisico, HornoIdeal
from motor.parametros import cargar_parametros
from motor.termopar import Termopar


def test_programa_rampa_y_meseta():
    p = cargar_parametros(BASE)  # 25 → 500 °C a 10 °C/min, meseta 30 min
    prog = Programa(p)
    assert prog.setpoint(0.0) == pytest.approx(298.15)
    assert prog.setpoint(60.0) == pytest.approx(308.15)
    assert prog.t_rampa == pytest.approx(47.5 * 60)
    assert prog.t_total == pytest.approx((47.5 + 30) * 60)
    assert prog.setpoint(prog.t_total) == pytest.approx(773.15)


def test_horno_ideal_exacto():
    h = HornoIdeal(200.0)
    T = 300.0
    for _ in range(50):
        T = h.paso(T, 700.0, 7.0)
    assert T == pytest.approx(700.0 + (300.0 - 700.0) * np.exp(-350.0 / 200.0), rel=1e-12)


def test_termopar_exacto():
    tc = Termopar(10.0)
    lectura = 300.0
    for _ in range(8):
        lectura = tc.paso(lectura, 500.0, 2.5)
    assert lectura == pytest.approx(500.0 - 200.0 * np.exp(-2.0), rel=1e-12)


def test_pid_saturacion():
    pid = PID(Kp=0.02, Ki=1e-3, Kd=0.0, t_muestreo=1.0)
    u = [pid.calcular(t, 800.0, 300.0) for t in np.arange(0.0, 200.0, 1.0)]
    assert u[-1] == 1.0
    I_antes = pid.integral
    pid.calcular(200.0, 800.0, 300.0)
    assert pid.integral == I_antes      # congelado mientras u satura
    # al bajar el error, sale de la saturación sin arrastre grande
    assert pid.calcular(201.0, 800.0, 799.0) < 1.0


def test_pid_muestreo():
    pid = PID(Kp=0.01, Ki=0.0, Kd=0.0, t_muestreo=5.0)
    u0 = pid.calcular(0.0, 400.0, 350.0)
    assert pid.calcular(2.0, 400.0, 300.0) == u0     # entre muestras mantiene la salida
    assert pid.calcular(5.0, 400.0, 300.0) == pytest.approx(1.0)


def test_ssr_onoff_promedio():
    t_ciclo, P, u = 2.0, 2500.0, 0.3
    for modo in ("onoff", "promedio"):
        ssr = SSR(modo, t_ciclo)
        E = sum(ssr.energia(t, t + 0.5, u, P) for t in np.arange(0.0, 10 * t_ciclo, 0.5))
        assert E == pytest.approx(u * P * 10 * t_ciclo, rel=1e-6), modo
    onoff = SSR("onoff", t_ciclo)
    assert onoff.potencia(0.1, u, P) == P and onoff.potencia(1.0, u, P) == 0.0
    assert SSR("promedio", t_ciclo).potencia(1.0, u, P) == pytest.approx(u * P)


def test_horno_fisico_estacionario():
    h = HornoFisico(C_h=2e4, UA_p=5.0, T_amb=298.15)
    T = 700.0
    P = 5.0 * (T - 298.15)
    assert h.paso(T, P, 0.0, 30.0) == pytest.approx(T, rel=1e-12)


def test_horno_fisico_balance():
    h = HornoFisico(C_h=2e4, UA_p=0.0, T_amb=298.15)
    assert h.paso(500.0, 1000.0, 200.0, 10.0) == pytest.approx(500.0 + 800.0 * 10.0 / 2e4)
