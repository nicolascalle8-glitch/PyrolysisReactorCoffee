"""Rutas comunes de las pruebas del motor."""
import os

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(RAIZ, "parametros", "base.toml")

import sys
# legacy/ y legacy/python al final de sys.path: solo las pruebas comparan contra legacy; el motor no lo importa
for _ruta in (os.path.join(RAIZ, "legacy"), os.path.join(RAIZ, "legacy", "python")):
    if _ruta not in sys.path:
        sys.path.append(_ruta)
