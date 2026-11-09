"""Rutas comunes de las pruebas del motor."""
import os

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(RAIZ, "parametros", "base.toml")

import sys
_LEGACY_PY = os.path.join(RAIZ, "legacy", "python")
if _LEGACY_PY not in sys.path:
    sys.path.append(_LEGACY_PY)  # al final: el motor nunca importa de legacy
