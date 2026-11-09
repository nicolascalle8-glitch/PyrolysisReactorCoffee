"""Hace importable el código legacy desde la raíz del repo: agrega legacy/ y legacy/python a sys.path."""
import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
for ruta in (AQUI, os.path.join(AQUI, "python")):
    if ruta not in sys.path:
        sys.path.insert(0, ruta)
