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
