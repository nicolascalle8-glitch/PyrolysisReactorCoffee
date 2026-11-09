"""
resultados.py — Lo que devuelve una simulación y cómo se exporta.

Resultado.rendimientos: % de la carga (masa húmeda inicial con cenizas): biochar (sólido carbonoso +
cenizas, sin la biomasa original que queda), bio_oil, gas, agua (humedad evaporada + agua de reacción),
sin_reaccionar (biomasa original que queda) y humedad_restante.
Resultado.historia: arreglos con t en min y T en °C. Resultado.mapas: lista de (t en min, matriz T r–z en °C).
balance_masa: error máximo del cierre de masa en el tiempo (% de la carga).
balance_energia: error acumulado máximo (% de la energía entregada al tubo).
"""

import csv
import json
import os
from dataclasses import dataclass, field

import numpy as np

TOL_MASA = 0.01      # % de la carga (spec §8)
TOL_ENERGIA = 0.5    # % de la energía entregada (spec §8)


@dataclass
class Resultado:
    rendimientos: dict
    historia: dict
    mapas: list
    balance_masa: float
    balance_energia: float
    supuestos: list = field(default_factory=list)
    r_mm: np.ndarray | None = None
    z_cm: np.ndarray | None = None
    n_lecho: int = 0
    avisos: list = field(default_factory=list)
    pasos: int = 0
    dT_max_paso: float = 0.0      # mayor ΔT aceptado en un paso (K)
    pasos_en_dt_min: int = 0      # pasos aceptados en el piso de dt sin cumplir los límites

    @property
    def confiable(self) -> bool:
        return self.balance_masa <= TOL_MASA and self.balance_energia <= TOL_ENERGIA

    def a_csv(self, carpeta: str) -> None:
        """Escribe historia.csv, rendimientos.json y mapas.npz en la carpeta (la crea si no existe)."""
        os.makedirs(carpeta, exist_ok=True)
        claves = list(self.historia)
        with open(os.path.join(carpeta, "historia.csv"), "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(claves)
            for fila in zip(*(self.historia[c] for c in claves)):
                w.writerow([f"{v:.6g}" for v in fila])
        with open(os.path.join(carpeta, "rendimientos.json"), "w", encoding="utf-8") as f:
            json.dump({"rendimientos": self.rendimientos, "balance_masa": self.balance_masa,
                       "balance_energia": self.balance_energia, "confiable": self.confiable,
                       "supuestos": self.supuestos, "avisos": self.avisos}, f, ensure_ascii=False, indent=2)
        np.savez(os.path.join(carpeta, "mapas.npz"),
                 t=np.array([t for t, _ in self.mapas]), T=np.array([m for _, m in self.mapas]),
                 r_mm=self.r_mm, z_cm=self.z_cm)
