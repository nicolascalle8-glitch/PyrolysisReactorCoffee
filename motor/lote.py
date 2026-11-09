"""
lote.py — Corre muchos casos en paralelo (estudios de sensibilidad, sintonía del PID).

simular_lote(lista, procesos=12) -> lista de Resultado o la excepción de cada caso, en el orden de entrada.
Usa ProcessPoolExecutor con un hilo BLAS por proceso (OMP/OPENBLAS/MKL/VECLIB_MAXIMUM_THREADS = 1) para
que los procesos no compitan por los núcleos. Un caso que falla no tumba el lote. procesos=1 corre en el
mismo proceso (útil para depurar).
"""

import os
from concurrent.futures import ProcessPoolExecutor

from .simulador import simular

VARIABLES_HILOS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")


def _correr(p):
    try:
        return simular(p)
    except Exception as e:  # el error queda en su posición del lote
        return e


def simular_lote(lista, procesos: int = 12):
    lista = list(lista)
    if procesos <= 1 or len(lista) <= 1:
        return [_correr(p) for p in lista]
    antes = {v: os.environ.get(v) for v in VARIABLES_HILOS}
    try:
        os.environ.update({v: "1" for v in VARIABLES_HILOS})  # lo heredan los procesos hijos
        with ProcessPoolExecutor(max_workers=min(procesos, len(lista))) as ex:
            return list(ex.map(_correr, lista))
    finally:
        for v, valor in antes.items():
            if valor is None:
                os.environ.pop(v, None)
            else:
                os.environ[v] = valor
