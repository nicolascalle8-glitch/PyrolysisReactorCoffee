# Spec — Motor nuevo del reactor vertical (subproyecto A)

Fecha: 2026-11-08 · Autor: Claude Opus 5.5 con Nick · Estado: **borrador para revisión de Nick**
Proyecto: horno-vertical (bitácora en Skill_Repository/proyectos/horno-vertical)

## 1. Para qué es

Nick quiere usar el modelo como un "tuner" de su horno real: calibrar el modelo con medidas del horno, sintonizar el PID y escoger la receta (rampa, T final, meseta). El código actual (`simulador_interactivo.py`) da resultados que cierran el balance, pero:

- mezcla menú, física, cinética vieja y gráficas en un solo archivo de ~600 líneas;
- tiene constantes escritas a mano sin fuente (malla, `h_comb = 90`, correlaciones de k y cp);
- modela el horno como un retardo que sigue el setpoint, sin PID, SSR ni potencia;
- usa un esquema explícito que el acero obliga a correr con pasos de ~0.3 s (25 s por corrida en un núcleo del M4).

Este spec cubre **solo el motor**. La calibración (B) y la app Streamlit (C) tendrán su propio spec y solo usarán la interfaz de la sección 7.

**Criterio de éxito:** el motor nuevo reproduce el caso base del modelo actual (sección 9), pasa las pruebas de verificación, corre al menos 5 veces más rápido, y cada parámetro tiene unidad, etiqueta y fuente.

## 2. Decisiones ya aprobadas (2026-11-08)

1. Orden del proyecto: A motor → B ensayo térmico y calibración térmica → C app Streamlit.
2. Numérica: paso dividido. Conducción implícita (Euler implícito, sistema disperso) + química integrada exacta por nodo.
3. Horno modelado como el real: masa térmica, resistencias, PID → SSR con ventana de tiempo, pérdidas, intercambio horno–tubo, termopar con retardo.
4. PID se sintoniza en Python (app); MATLAB/Simulink verifica y sirve para el microcontrolador.
5. Repo PyrolysisReactorCoffee ordenado como en la sección 3; código actual a `legacy/`; datos del GIEM y FTIR **no** se suben (repo público, el reporte trae datos personales).

## 3. Estructura del repo

```
PyrolysisReactorCoffee/
├── README.md            Qué es, cómo instalar, cómo correr (1 página, rutas de Mac)
├── pyproject.toml       Dependencias; reemplaza requirements.txt (está en UTF-16 y pip no lo lee bien)
├── parametros/
│   └── base.toml        Caso base con unidad, etiqueta, fuente y rango por parámetro
├── motor/               Paquete del motor (sección 4)
├── calibracion/         Subproyecto B (vacío con README de una línea)
├── app/                 Subproyecto C (vacío con README de una línea)
├── matlab/              Verificación y planta para Simulink (se mueve lo que hay)
├── tests/               Pruebas del motor
├── docs/                Specs, EVALUACION_CINETICAS.md, figuras
└── legacy/              Código actual sin cambios de lógica, movido con git mv
```

- `legacy/` contiene: `simulador_interactivo.py`, `python/`, `run_simulation.py`, `run_lumped.py`, `tests/` actuales. Solo se ajustan rutas de import para que sus pruebas sigan corriendo desde `legacy/`.
- `legacy/` se borra solo cuando el motor pase la sección 9 y Nick lo apruebe. **Hecho el 2026-11-08** (aprobado por Nick); queda en el historial de git (commit c2acaae).
- `.gitignore` incluye `resultados/` y `datos/` para que nada de laboratorio se suba por error.
- Python ≥ 3.10. TOML con `tomllib` (3.11+) o `tomli` (3.10).

## 4. Módulos del motor

| Módulo | Qué hace | Depende de |
|---|---|---|
| `parametros.py` | Lee `base.toml`, aplica cambios, valida unidades, rangos y composición; lista los [Supuesto] | — |
| `malla.py` | Nodos r–z del lecho y del acero, volúmenes de control, áreas de cara | parametros |
| `propiedades.py` | k, ρ, cp efectivas del lecho según T y conversión; acero constante | parametros |
| `horno.py` | Balance de energía del horno (modo físico) o retardo de primer orden (modo ideal, para comparar con legacy) | parametros |
| `control.py` | Programa rampa/meseta, PID discreto con anti-windup, SSR por ventana | — |
| `termopar.py` | Lectura con retardo de primer orden en una posición (r, z) | malla |
| `conduccion.py` | Arma y resuelve A·Tⁿ⁺¹ = b (disperso) para lecho + acero | malla, propiedades |
| `cinetica/rac.py` | Esquema RAC + TGL (se reutiliza `kinetics_rac.py` sin cambiar la química) | — |
| `cinetica/secado.py` | Evaporación de la humedad | parametros |
| `simulador.py` | Bucle del paso (sección 5); sin física adentro | todos |
| `resultados.py` | Historias, rendimientos, balances, exportar CSV/JSON | — |
| `graficas.py` | Gráficas (temperaturas, masa, mapas 2D); la app las reutiliza | resultados |

Cada módulo tiene un docstring corto que dice qué hace, qué recibe, qué devuelve y de dónde salen sus ecuaciones.

## 5. Un paso de tiempo

En cada paso `dt`:

1. **Control.** Setpoint del programa. El PID usa la lectura del termopar de control y entrega el ciclo útil u ∈ [0, 1].
   - PID en forma estándar: u = Kp·e + Ki·∫e dt + Kd·de/dt, con muestreo `t_muestreo` y anti-windup por saturación. La conversión a las unidades del MC-5438 (banda P, etc.) se agrega cuando llegue el manual [Por validar].
   - Modo `promedio`: potencia = P_max·u en todo el paso.
   - Modo `onoff`: el paso no pasa del tiempo de ciclo del SSR; dentro de cada ventana las resistencias están prendidas u·t_ciclo y apagadas el resto.
2. **Horno (modo físico).**
   C_h·dT_h/dt = P − UA_p·(T_h − T_amb) − Q_tubo
   Q_tubo = Σ_z A_z·[h_c·(T_h − T_w,z) + ε·σ·(T_h⁴ − T_w,z⁴)]
   Se integra con sub-pasos explícitos si hace falta (es una sola ecuación).
   **Modo ideal:** T_h sigue al setpoint con un retardo τ_piedra (igual que legacy).
3. **Calor en el reactor (implícito).** Euler implícito sobre lecho + acero:
   - Interior en coordenadas cilíndricas con volúmenes finitos; en r = 0 el volumen central (equivale a L'Hôpital).
   - Interfaz lecho–acero: resistencia de contacto en serie dentro de la matriz → el flujo que sale del acero es el que entra al lecho.
   - Cara externa del acero: Robin con h_c + radiación linealizada, h_rad = ε·σ·(T_h² + T_w²)·(T_h + T_w).
   - Caras superior e inferior: adiabáticas [Supuesto, igual que legacy]; se puede cambiar por parámetro.
   - k y cp evaluados con la T del paso anterior. Solver: `scipy.sparse.linalg.splu` o `spsolve`.
4. **Química.** Con Tⁿ⁺¹, cada nodo reacciona con la integración exacta actual (B·(1 − e^(−k·dt))). Calor de reacción y de evaporación → ajuste de temperatura ΔT = q·dt/(ρ·cp). Masa y energía cierran por construcción.
5. **Control de `dt`.** Inicial `numerica.dt0` (2 s [Supuesto]). Si algún nodo cambia más de `numerica.dT_max` (5 K [Supuesto]) o la conversión de un nodo cambia más de `numerica.dX_max` (5 % [Supuesto]), se repite el paso con dt/2. Si el cambio máximo es menor que la mitad del límite, dt crece ×1.5 hasta `numerica.dt_max`. En modo `onoff`, dt ≤ t_ciclo.
6. **Registro** cada `numerica.t_registro` s: temperaturas en nodos clave y termopares, masas, balances.

## 6. Parámetros (`parametros/base.toml`)

Cada parámetro:

```toml
[horno.P_max]
valor    = 2500
unidad   = "W"
etiqueta = "Supuesto"          # Dato | Medido | Norma | Supuesto | Calculado | Por validar
fuente   = "config.py original, sin fuente"
min      = 500
max      = 5000
```

Grupos: `programa`, `control`, `horno`, `geometria`, `malla`, `lecho`, `cinetica`, `termopares`, `numerica`.

El caso base copia los valores de legacy para poder compararlos (todos [Supuesto] salvo los que tienen fuente):

| Parámetro | Valor base | Etiqueta / fuente |
|---|---|---|
| Rampa, T final, meseta | 10 °C/min, 500 °C, 30 min | Caso base del informe EVALUACION_CINETICAS |
| R lecho, espesor acero, L | 25 mm, 2.5 mm, 0.20 m | Supuesto (legacy) |
| Nr lecho, Nr acero, Nz | 14, 4, 20 | Supuesto (legacy) → estudio de malla |
| Densidad aparente | 380 kg/m³ | Supuesto (legacy) |
| Humedad de carga | 8 % bh | Supuesto (legacy; el GIEM midió 61.2 % bh en fresco) |
| Cenizas | 2.31 % bs | **Dato**: GIEM 28ago2604 |
| Celulosa / hemicelulosa / lignina / lípidos | 12 / 35 / 28 / 15 % bs | **Supuesto sin fuente** (no está en el reporte GIEM) |
| k, cp del lecho | correlaciones lineales de legacy | Supuesto (legacy) |
| h horno–tubo | 90 W/m²K | Supuesto (legacy `h_comb`) |
| Resistencia extra de contacto | 0.003 m²K/W | Supuesto (legacy) |
| τ piedra (modo ideal) | 200 s | Supuesto (legacy) |
| P_max | 2500 W | Supuesto (config.py) |
| C_h, UA_p, ε | por calibrar | Por validar (subproyecto B) |
| Kp, Ki, Kd, t_ciclo | por definir | Por validar (manual MC-5438) |

## 7. Interfaz pública (lo único que usan B, C y Cursor)

```python
from motor import cargar_parametros, simular, simular_lote

p = cargar_parametros("parametros/base.toml", cambios={"control.Kp": 4.0})
r = simular(p)
r.rendimientos      # dict: biochar, bio_oil, gas, agua, sin_reaccionar, humedad_restante (% de la carga)
r.historia          # dict de arrays: t, T_horno, T_set, T_pared, T_centro, T_termopar_<n>, masa_residual
r.mapas             # lista de (t, matriz T r–z) en instantes clave
r.balance_masa      # error máximo (% de la carga)
r.balance_energia   # error máximo (% de la energía entregada)
r.confiable         # False si algún balance pasa la tolerancia
r.a_csv("resultados/<nombre>/")

lote = simular_lote([p1, p2, ...], procesos=12)   # ProcessPoolExecutor, 1 hilo BLAS por proceso
```

## 8. Errores

- Unidad distinta a la esperada, valor fuera de [min, max], composición que no suma 100 % (con "otros"): `ErrorParametro` con el nombre del parámetro y el motivo. No corre.
- Al cargar: lista de parámetros [Supuesto] y [Por validar] (la app la muestra).
- T no finita o divergente: `ErrorSimulacion` con el tiempo y el nodo.
- Balance de masa > 0.01 % o de energía > 0.5 %: el resultado sale con `confiable = False` y un aviso; no se oculta.
- En `simular_lote`, un caso que falla no tumba el lote: queda registrado con su error.

## 9. Verificación (segundo camino)

| # | Prueba | Criterio de aceptación |
|---|---|---|
| 1 | Cilindro sin reacción, pared a T fija, contra la solución exacta en serie de Bessel (conducción transitoria en cilindro; Çengel/Incropera, edición a confirmar con el libro de Nick) | Error máximo ≤ 1 % de la diferencia de temperatura inicial |
| 2 | Caso base en modo horno ideal contra legacy (biochar 24.93, bio-oil 38.94, gas 20.23, agua 15.90 %) | Cada rendimiento dentro de ±1.0 punto porcentual; curva T_centro con error RMS ≤ 10 K [Criterio, se ajusta si hay razón física] |
| 3 | Balance de masa y energía en los 16 casos de las pruebas actuales | Masa ≤ 0.01 %; energía ≤ 0.5 % |
| 4 | Convergencia en dt: dt0 y dt0/2 | Rendimientos cambian < 0.2 puntos |
| 5 | Convergencia de malla: malla base y malla ×2 en r y z | Rendimientos cambian < 0.5 puntos; se reporta el rango |
| 6 | Velocidad: caso base en un núcleo del M4 | ≥ 5 veces más rápido que legacy (25 s) |
| 7 | Parte térmica en MATLAB (ode15s sobre la misma malla, sin reacción) | Mismo criterio que la prueba 1 entre Python y MATLAB |
| 8 | Validación de parámetros: casos con unidad mala, fuera de rango y composición que no suma | Cada uno lanza `ErrorParametro` con el nombre correcto |

Si una prueba no pasa, no se entrega hasta saber por qué.

## 10. Fuera de este spec

- Calibración contra termopares (B) y la app Streamlit (C).
- Craqueo de alquitranes en fase gas (se decide cuando haya rendimientos medidos).
- Cinética de proteínas (el esquema RAC no la tiene; la proteína entra repartida en "otros" hasta tener fuente).
- Sintonía en el microcontrolador real.

## 11. Riesgos

- **Identificabilidad** (lo dijo el panel): con solo termopares, h de contacto y k del lecho se confunden. Se mitiga en B con dos termopares (pared y centro) y una corrida en vacío.
- **Error del paso dividido en la exoterma** (min 40–45): se controla con el dt adaptativo y la prueba 4.
- **Prueba 2 puede no pasar exacta** porque legacy usa otra numérica. Si la diferencia viene del esquema explícito de legacy (paso más grueso o error de interfaz), se documenta con la prueba 1 y la 4 como árbitros y se le muestra a Nick antes de cambiar el criterio.
