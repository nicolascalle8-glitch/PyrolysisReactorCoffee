# Reactor vertical de pirólisis de residuo de café

Modelo 2D (r, z) del reactor vertical del horno (práctica GIEM, UdeA/EIA): lecho de café dentro de un tubo
de acero inoxidable, calentado por un horno con resistencias, PID y SSR. Cinética RAC + triglicéridos.
El motor (`motor/`) es la base del "tuner" del horno: calibrar el modelo, sintonizar el PID y escoger la receta.

## Instalar (Mac)

```bash
cd ~/Documents/Documents\ -\ MacBook\ Pro\ de\ Nicolas/Max_Power/01_Horno_Vertical/modelos/PyrolysisReactorCoffee
python3 -m venv .venv
.venv/bin/pip install -e ".[pruebas]"
```

## Carpetas

| Carpeta | Qué hay |
|---|---|
| `motor/` | Paquete del motor nuevo (parámetros, malla, conducción implícita, cinética, horno, simulador) |
| `parametros/base.toml` | Caso base: cada parámetro con unidad, etiqueta, fuente y rango |
| `tests/` | Pruebas del motor (`python -m pytest tests -q`) |
| `calibracion/` | Subproyecto B (vacío) |
| `app/` | Subproyecto C, app Streamlit (vacío) |
| `matlab/` | Modelos MATLAB y verificación térmica |
| `docs/` | Specs, planes, evaluación de cinéticas |
| `legacy/` | Código anterior, sin cambios de lógica (se borra con aprobación de Nick) |

`datos/` y `resultados/` están en `.gitignore`: nada de laboratorio se sube (el repo es público).

## Correr un caso

```python
from motor import cargar_parametros, simular

p = cargar_parametros("parametros/base.toml", cambios={"programa.rampa": 5.0})
r = simular(p)
print(r.rendimientos)            # % de la carga: biochar, bio_oil, gas, agua, sin_reaccionar, humedad_restante
print(r.balance_masa, r.balance_energia, r.confiable)
r.a_csv("resultados/rampa_5/")   # historia.csv, rendimientos.json, mapas.npz

from motor.graficas import graficar
graficar(r, "resultados/rampa_5/")
```

`p.supuestos()` lista los parámetros [Supuesto] y [Por validar]. Una clave mal escrita, una unidad distinta o un
valor fuera de rango lanzan `ErrorParametro` con el nombre del parámetro.

## Correr un lote

```python
from motor import cargar_parametros, simular_lote

base = cargar_parametros("parametros/base.toml")
casos = [base.con_cambios({"programa.rampa": v}) for v in (2, 5, 10, 20, 50)]
resultados = simular_lote(casos, procesos=12)   # un caso que falla queda como excepción en su posición
```

En macOS el lote necesita el `if __name__ == "__main__":` en el script (los procesos se crean con *spawn*).

## Pruebas

```bash
.venv/bin/python -m pytest tests -q              # todas (≈3 min)
.venv/bin/python -m pytest tests -q -m "not lento"  # rápidas
```

Las pruebas 1–8 del spec (`docs/specs/2026-11-08-motor-nuevo-design.md` §9) están en `tests/`. La prueba 7
(MATLAB) está en `matlab/verificacion/`.
