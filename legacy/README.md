# legacy/

Código anterior al motor nuevo, movido sin cambios de lógica (2026-11-08). Sirve de referencia para la
prueba 2 (caso base) y se borra solo cuando Nick lo apruebe.

- `simulador_interactivo.py`: modelo vigente hasta el motor nuevo (explícito, horno como retardo).
- `python/`: modelos anteriores y la cinética RAC + TGL (`kinetics_rac.py`).
- Pruebas: `python -m pytest legacy/tests -q` (desde la raíz del repo).
