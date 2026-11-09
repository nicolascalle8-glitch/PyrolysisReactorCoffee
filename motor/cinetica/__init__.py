"""Cinética del lecho: esquema RAC + TGL (rac.py) y secado (secado.py)."""
from .rac import FRACCIONES, REGIMENES, CineticaRAC
from .secado import paso_secado

__all__ = ["CineticaRAC", "REGIMENES", "FRACCIONES", "paso_secado"]
