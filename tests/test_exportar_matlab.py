"""Prueba del exportador para la verificación térmica en MATLAB (prueba 7 del spec)."""
import csv
import os
import sys

import numpy as np

from conftest import RAIZ

sys.path.insert(0, os.path.join(RAIZ, "matlab", "verificacion"))
import exportar_caso  # noqa: E402


def test_exporta_csv(tmp_path):
    rutas = exportar_caso.exportar(str(tmp_path))
    nodos, meta = rutas
    assert os.path.exists(nodos) and os.path.exists(meta)
    with open(nodos, encoding="utf-8") as f:
        cab = next(csv.reader(f))
    assert cab[:8] == ["i", "k", "r", "z", "V", "k_W_mK", "rho_cp", "T0"]
    datos = np.loadtxt(nodos, delimiter=",", skiprows=1)
    assert datos.shape[0] == 18 * 20 and np.all(np.isfinite(datos))
    with open(meta, encoding="utf-8") as f:
        m = {fila[0]: fila[1] for fila in csv.reader(f) if fila and fila[0] != "nombre"}
    assert float(m["T_horno"]) > float(m["T0"])
    tiempos = [float(c[2:]) for c in cab if c.startswith("T_")]
    assert len(tiempos) >= 10 and tiempos[0] > 0


def test_segundo_camino_bdf(tmp_path):
    """Euler implícito del motor contra el sistema semidiscreto integrado con BDF (lo mismo que ode15s)."""
    sol = exportar_caso.verificar_con_bdf(str(tmp_path))
    assert sol["error_max_K"] <= 0.01 * sol["dT0_K"], sol
