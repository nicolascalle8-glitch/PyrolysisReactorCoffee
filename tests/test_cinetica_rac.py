"""
test_cinetica_rac.py
Pruebas de la cinética RAC + triglicéridos (python/kinetics_rac.py):
  - cada reacción transcrita está balanceada en C, H y O (detecta errores al copiar coeficientes),
  - la masa se conserva para cualquier paso de tiempo,
  - reproduce los datos experimentales de café usados como validación.

Ejecutar desde la raíz del proyecto:
    python -m pytest tests -v
"""

import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python"))

import comparar_cineticas as cc
from kinetics_rac import CineticaRAC, REGIMENES


@pytest.mark.parametrize("regimen", REGIMENES)
def test_reacciones_balanceadas_atomo_a_atomo(regimen):
    kin = CineticaRAC(regimen)
    desbalanceadas = {n: e for n, e in kin.error_atomico.items() if e > 1e-6}
    assert not desbalanceadas, f"Reacciones con desbalance C/H/O: {desbalanceadas}"
    assert np.allclose(kin.Y.sum(axis=0), 1.0)


@pytest.mark.parametrize("dt", [1e-3, 1.0, 1e3])
@pytest.mark.parametrize("T", [500.0, 700.0, 1100.0])
def test_masa_se_conserva_para_cualquier_paso(T, dt):
    kin = CineticaRAC('lecho_fijo')
    m = kin.estado_inicial(1.0, forma=(3, 4))
    total = {'gas': 0.0, 'agua_rxn': 0.0, 'bio_oil': 0.0}
    for _ in range(20):
        _, lib = kin.paso(m, np.full((3, 4), T), dt)
        for f in total:
            total[f] += lib[f]
    assert np.all(m >= 0.0)
    balance = kin.masa_solida(m) + sum(total.values())
    assert np.allclose(balance, 1.0, atol=1e-12)


def _rac(caso):
    T_de_t, t_total, _ = cc.programa(caso['beta'], caso['T_final'], caso['meseta'])
    return cc.simular_rac(caso['regimen'], T_de_t, t_total, cc.dt_programa(caso['beta']), caso['w_H2O'])


def test_validacion_pirolisis_rapida_bok_2012():
    """Bok et al. (2012): máximo de bio-oil 54.85 % a 823 K."""
    r = _rac(cc.CASOS[1])
    assert abs(100 * (r['bio_oil'] + r['agua']) - 54.85) < 3.0


def test_validacion_lecho_fijo_inayat_2022_liquido():
    """Inayat et al. (2022): bio-oil ~22.5 % a 300 °C, 30 min, ~0.9 °C/min."""
    r = _rac(cc.CASOS[0])
    assert abs(100 * (r['bio_oil'] + r['agua']) - 22.5) < 5.0


def test_rapida_produce_mas_liquido_que_lenta():
    rapida = _rac(cc.CASOS[1])
    lenta = _rac(dict(cc.CASOS[1], beta=10.0, meseta=30.0, regimen='lecho_fijo'))
    assert rapida['bio_oil'] > lenta['bio_oil']
    assert rapida['biochar'] < lenta['biochar']
