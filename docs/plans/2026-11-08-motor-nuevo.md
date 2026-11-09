# Motor nuevo — Plan de implementación

> **Para agentes:** SUB-SKILL REQUERIDA: superpowers:executing-plans (recomendado aquí) o superpowers:subagent-driven-development. Pasos con casillas (`- [ ]`).

**Goal:** Reemplazar `simulador_interactivo.py` por un paquete `motor/` modular, implícito en la conducción, con horno físico + PID/SSR y parámetros con fuente, que reproduzca el caso base y corra ≥ 5 veces más rápido.

**Architecture:** Paso dividido: control → horno → conducción implícita (matriz dispersa) → química exacta por nodo, con dt adaptativo. Parámetros en `parametros/base.toml`. Interfaz pública: `cargar_parametros`, `simular`, `simular_lote`.

**Tech Stack:** Python ≥ 3.10, numpy, scipy (sparse), matplotlib, tomllib/tomli, pytest. MATLAB solo para la prueba 7.

**Spec:** `docs/specs/2026-11-08-motor-nuevo-design.md` (aprobado por Nick 2026-11-08).

## Global Constraints

- Python ≥ 3.10; dependencias solo numpy, scipy, matplotlib, tomli (si < 3.11); pytest para pruebas.
- Nombres, docstrings y mensajes de error en español simple.
- `legacy/` no cambia de lógica; solo rutas de import.
- Nada de datos de laboratorio en el repo (`datos/` y `resultados/` en `.gitignore`).
- Tolerancias del spec §9: Bessel ≤ 1 % de ΔT₀; caso base ±1.0 pp por rendimiento y RMS T_centro ≤ 10 K; masa ≤ 0.01 %; energía ≤ 0.5 %; dt/2 cambia < 0.2 pp; malla ×2 cambia < 0.5 pp; velocidad ≥ 5× legacy (25 s → ≤ 5 s en un núcleo M4).
- Valores numéricos de partida: dt0 = 2 s, dT_max = 5 K, dX_max = 0.05, crecimiento ×1.5, reducción ×0.5.
- Comandos desde la raíz del repo: `python -m pytest tests -q`.

## Review Focus

1. **T final por debajo de la pirólisis (p. ej. 150 °C):** corre, conversión ~0, balances cierran, sin divisiones por cero. → Tarea 8.
2. **Humedad 0 %:** el secado no divide por cero y el resultado es válido. → Tarea 6.
3. **Clave inexistente en `cambios` (error de tipeo desde la app):** `ErrorParametro` que nombra la clave. → Tarea 2.
4. **Rampa muy rápida (50 °C/min) y meseta 0 min:** el dt adaptativo no se queda en bucle y termina. → Tarea 9.
5. **Malla muy gruesa (Nr lecho = 3) o fuera de rango:** corre si está en [min, max]; si no, `ErrorParametro`. → Tareas 2 y 3.

---

### Task 1: Ordenar el repo

**Files:**
- Move (git mv): `simulador_interactivo.py`, `python/`, `run_simulation.py`, `run_lumped.py`, `tests/` → `legacy/`
- Create: `pyproject.toml`, `.gitignore`, `README.md` (reescrito), `calibracion/README.md`, `app/README.md`, `legacy/README.md`, `legacy/conftest.py`
- Delete: `requirements.txt` (UTF-16), `pyrightconfig.json` si apunta a rutas viejas
- Move: `docs/EVALUACION_CINETICAS.*` y `docs/figuras/` se quedan en `docs/`

**Interfaces:** Produces: estructura del spec §3; `legacy/` importable con `legacy/conftest.py` que agrega `legacy/` y `legacy/python` a `sys.path`.

- [ ] **Step 1:** `git mv` de los archivos a `legacy/`.
- [ ] **Step 2:** Crear `legacy/conftest.py` (ajuste de `sys.path`) y corregir en `legacy/tests/*.py` la variable `ROOT` para que apunte a `legacy/`.
- [ ] **Step 3:** Run `python -m pytest legacy/tests -q -k "base or rac"` → PASS (prueba rápida de que legacy sigue vivo).
- [ ] **Step 4:** `pyproject.toml` (proyecto `pirolisis-cafe`, paquete `motor`, deps del Global Constraints, `[tool.pytest.ini_options] testpaths=["tests"]`), `.gitignore` (`resultados/`, `datos/`, `.venv/`, `__pycache__/`), README de 1 página con instalación en Mac (`python3 -m venv .venv && .venv/bin/pip install -e .`) y mapa de carpetas; READMEs de una línea en `calibracion/` y `app/`.
- [ ] **Step 5:** Commit `refactor: ordena el repo (legacy/, pyproject, README)`.

### Task 2: Parámetros

**Files:** Create `parametros/base.toml`, `motor/__init__.py`, `motor/parametros.py`, `tests/test_parametros.py`

**Interfaces:**
- Produces: `class ErrorParametro(ValueError)`; `@dataclass Parametro(valor: float|str, unidad: str, etiqueta: str, fuente: str, min: float|None, max: float|None)`; `class Parametros` con `p["grupo.nombre"] -> float|str` (valor), `p.info("grupo.nombre") -> Parametro`, `p.supuestos() -> list[str]` (claves con etiqueta Supuesto o Por validar), `p.con_cambios(dict) -> Parametros`; `cargar_parametros(ruta: str, cambios: dict|None = None) -> Parametros`.
- Grupos y claves exactas (todas en base.toml con los valores de la tabla del spec §6): `programa.{T_inicial, rampa, T_final, meseta}` (°C, °C/min, °C, min); `control.{modo_ssr, Kp, Ki, Kd, t_ciclo, t_muestreo}`; `horno.{modo, tau_piedra, P_max, C_h, UA_p, h_c, emisividad, T_amb}`; `geometria.{R_lecho, e_acero, L}`; `malla.{Nr_lecho, Nr_acero, Nz}`; `lecho.{rho_aparente, humedad, cenizas, celulosa, hemicelulosa, lignina, lipidos, otros, R_contacto_extra}`; `acero.{rho, cp, k}`; `cinetica.{regimen, acetilos_hce}`; `termopares.{control_r, control_z, tau}`; `numerica.{dt0, dt_max, dT_max, dX_max, t_registro}`.
- `horno.modo = "ideal"` y `control.modo_ssr = "promedio"` en base.toml (para la prueba 2). Composición: `otros = 100 − suma` de las demás orgánicas en base seca, de modo que celulosa+hemicelulosa+lignina+lipidos+otros+cenizas = 100.

- [ ] **Step 1: Tests** en `tests/test_parametros.py`:
  - `test_carga_base`: `p["programa.rampa"] == 10.0`, `p.info("lecho.cenizas").etiqueta == "Dato"`, `"lecho.celulosa" in p.supuestos()`.
  - `test_cambio_valido`: `cargar_parametros(B, {"control.Kp": 4.0})["control.Kp"] == 4.0`.
  - `test_clave_inexistente`: `{"control.kp": 4}` → `ErrorParametro` cuyo mensaje contiene `"control.kp"`.
  - `test_fuera_de_rango`: `{"malla.Nr_lecho": 1}` → `ErrorParametro` con `"malla.Nr_lecho"`; `{"malla.Nr_lecho": 3}` carga bien.
  - `test_unidad_mala`: TOML temporal con `programa.rampa` en `"K/s"` → `ErrorParametro` con `"unidad"`.
  - `test_composicion_no_suma`: `{"lecho.celulosa": 50}` → `ErrorParametro` con `"composición"`.
- [ ] **Step 2:** Run `python -m pytest tests/test_parametros.py -q` → FAIL (módulo no existe).
- [ ] **Step 3:** Implementar `motor/parametros.py` y `base.toml`. Unidades esperadas en un dict `UNIDADES` dentro de `parametros.py`; tolerancia de suma de composición 0.01.
- [ ] **Step 4:** Run → PASS.
- [ ] **Step 5:** Commit `feat(motor): parámetros con unidad, etiqueta, fuente y validación`.

### Task 3: Malla

**Files:** Create `motor/malla.py`, `tests/test_malla.py`

**Interfaces:** Consumes `Parametros`. Produces `@dataclass Malla` con `r: ndarray (Nr_total,)`, `z: ndarray (Nz,)`, `n_lecho: int`, `es_acero: ndarray[bool]`, `V: ndarray (Nr_total, Nz)` volúmenes de control (m³), `A_r: ndarray (Nr_total+1,)` áreas de cara radial por unidad de dz (2π r_cara), `dz: float`, `R_lecho`, `R_ext`; `crear_malla(p: Parametros) -> Malla`. Nodo 0 en r = 0 (volumen central π(dr/2)²·dz), igual que legacy.

- [ ] **Step 1: Tests:** `test_volumen_total`: `V[lecho].sum() ≈ π R_lecho² L` y `V[acero].sum() ≈ π (R_ext² − R_lecho²) L` (rel 1e-12); `test_malla_gruesa`: Nr_lecho = 3 crea malla válida.
- [ ] **Step 2:** Run → FAIL. **Step 3:** Implementar. **Step 4:** Run → PASS.
- [ ] **Step 5:** Commit `feat(motor): malla r–z de lecho y acero`.

### Task 4: Propiedades

**Files:** Create `motor/propiedades.py`, `tests/test_propiedades.py`

**Interfaces:** Produces `k_lecho(T: ndarray, conv: ndarray) -> ndarray`, `cp_lecho(T, conv) -> ndarray`, `rho_lecho(masas: ndarray) -> ndarray` (piso 50 kg/m³), con las correlaciones de legacy (spec §6, `simulador_interactivo.py` líneas de `k_bed`, `cp_bed`, `rho_bed`), T en K.

- [ ] **Step 1: Test** `test_igual_a_legacy`: para T ∈ {300, 500, 800} K y conv ∈ {0, 0.5, 1}, valores iguales a las fórmulas de legacy (rel 1e-12).
- [ ] **Step 2–4:** FAIL → implementar → PASS. **Step 5:** Commit.

### Task 5: Conducción implícita

**Files:** Create `motor/conduccion.py`, `tests/test_conduccion.py`

**Interfaces:** Consumes `Malla`, propiedades. Produces `paso_conduccion(malla, T: ndarray, k: ndarray, rho_cp: ndarray, h_ext: ndarray (Nz,), T_horno: float, R_contacto: ndarray (Nz,), dt: float, adiabatico_z: bool = True) -> tuple[ndarray, float]` → (Tⁿ⁺¹, calor que entró por la cara externa en J). Euler implícito, volúmenes finitos, `scipy.sparse.csr_matrix` + `scipy.sparse.linalg.spsolve`. Resistencia lecho–acero en serie: `dr_l/(2k_l) + dr_a/(2k_a) + R_contacto_extra`. Radiación entra por `h_ext` (la calcula el simulador).

- [ ] **Step 1: Tests:**
  - `test_bessel` (prueba 1): solo lecho homogéneo (k, ρcp constantes), h_ext = 1e9, T₀ = 300 K, T_horno = 400 K, a Fo = 0.05, 0.2, 0.5 comparar T(r) contra la serie θ = Σ 2/(λₙ J₁(λₙ)) J₀(λₙ r/R) e^(−λₙ² Fo) con 50 términos (`scipy.special.jn_zeros`, `j0`, `j1`); error máx ≤ 0.01·(400−300).
  - `test_energia_adiabatica`: h_ext = 0 → Σ ρcp·V·T se conserva (rel 1e-10).
  - `test_flujo_interfaz`: con h_ext > 0, calor entrado = Δ(Σ ρcp·V·T) (rel 1e-9).
- [ ] **Step 2:** FAIL. **Step 3:** Implementar. **Step 4:** PASS.
- [ ] **Step 5:** Commit `feat(motor): conducción implícita con verificación de Bessel`.

### Task 6: Cinética y secado

**Files:** Create `motor/cinetica/__init__.py`, `motor/cinetica/rac.py` (copia de `legacy/python/kinetics_rac.py` sin cambiar la química), `motor/cinetica/secado.py`, `tests/test_cinetica.py`

**Interfaces:** Produces `CineticaRAC` (misma API que legacy: `estado_inicial`, `paso(m, T, dt) -> (q, liberado)`, `masa_solida`, `masa_original`); `paso_secado(agua: ndarray, T: ndarray, dt: float) -> tuple[agua_nueva, evaporada, q]` con A = 5.0e5 s⁻¹, Ea = 41.8e3 J/mol, ΔH = −2.26e6 J/kg y umbral T > 350 K (valores de legacy, [Supuesto]).

- [ ] **Step 1: Tests:** traer `legacy/tests/test_cinetica_rac.py` adaptado a `motor.cinetica` (balance atómico, conservación); `test_secado_humedad_cero`: agua = 0 → sin NaN, evaporada = 0; `test_secado_conserva`: agua_nueva + evaporada = agua (rel 1e-12).
- [ ] **Step 2–4:** FAIL → implementar → PASS. **Step 5:** Commit.

### Task 7: Horno, control y termopar

**Files:** Create `motor/horno.py`, `motor/control.py`, `motor/termopar.py`, `tests/test_horno_control.py`

**Interfaces:** Produces:
- `Programa(p).setpoint(t: float) -> float` (K).
- `PID(Kp, Ki, Kd, t_muestreo)` con `.calcular(t, setpoint, medida) -> float` en [0, 1], anti-windup por saturación (integrador congelado cuando u satura).
- `SSR(modo: str, t_ciclo: float)` con `.potencia(t, u, P_max) -> float` (promedio: u·P_max; onoff: P_max si `(t mod t_ciclo) < u·t_ciclo`, si no 0).
- `HornoIdeal(tau)` con `.paso(T_h, T_set, dt) -> T_h` (exacto: T_set + (T_h − T_set)·e^(−dt/τ)); `HornoFisico(C_h, UA_p, T_amb)` con `.paso(T_h, P, Q_tubo, dt) -> T_h`.
- `Termopar(tau)` con `.paso(T_lectura, T_real, dt) -> T_lectura` (exacto, primer orden).

- [ ] **Step 1: Tests:** horno ideal contra la solución exacta; termopar igual; `test_pid_saturacion`: error grande sostenido → u = 1 y el integrador no crece más; `test_ssr_onoff_promedio`: en 10 ciclos con u = 0.3 la energía entregada = 0.3·P_max·10·t_ciclo (rel 1e-6); `test_horno_fisico_estacionario`: P = UA_p·(T − T_amb), Q_tubo = 0 → T constante.
- [ ] **Step 2–4:** FAIL → implementar → PASS. **Step 5:** Commit.

### Task 8: Simulador y resultados

**Files:** Create `motor/simulador.py`, `motor/resultados.py`, `tests/test_simulador.py`; Modify `motor/__init__.py` (exporta `cargar_parametros`, `simular`, `ErrorParametro`, `ErrorSimulacion`)

**Interfaces:** Consumes todo lo anterior. Produces `class ErrorSimulacion(RuntimeError)`; `simular(p: Parametros) -> Resultado`; `@dataclass Resultado` con `rendimientos: dict` (claves `biochar, bio_oil, gas, agua, sin_reaccionar, humedad_restante`, en %), `historia: dict[str, ndarray]` (claves `t, T_horno, T_set, T_pared, T_centro, T_termopar_control, masa_residual`; T en °C, t en min), `mapas: list[tuple[float, ndarray]]`, `balance_masa: float`, `balance_energia: float`, `confiable: bool`, `supuestos: list[str]`, `.a_csv(carpeta: str) -> None` (historia.csv, rendimientos.json, mapas.npz). Orden del paso: spec §5. Primero con dt fijo = dt0 (el adaptativo va en la Tarea 9).

- [ ] **Step 1: Tests:**
  - `test_caso_base_vs_legacy` (prueba 2): base.toml → cada rendimiento a ±1.0 pp de {biochar 24.93, bio_oil 38.94, gas 20.23, agua 15.90}; RMS de T_centro contra legacy (correr `legacy.simulador_interactivo.ejecutar_simulacion` con el cfg equivalente, interpolado al mismo t) ≤ 10 K. Marcar `@pytest.mark.lento`.
  - `test_balances_16_casos` (prueba 3): los 12 casos de `CASOS_INTERACTIVO` de legacy traducidos a `cambios` + 4 de `simulation.py` (rampa 5 y 20, R 60 mm, humedad 40 %) → `balance_masa ≤ 0.01` y `balance_energia ≤ 0.5`.
  - `test_temperatura_baja`: T_final = 150 °C → corre, `rendimientos["sin_reaccionar"] > 90`, `confiable`.
  - `test_divergencia`: monkeypatch de `paso_conduccion` para que devuelva un `NaN` en un nodo → `ErrorSimulacion` cuyo mensaje contiene `"t ="` y `"nodo"`.
  - `test_a_csv`: escribe los 3 archivos en `tmp_path`.
- [ ] **Step 2:** FAIL. **Step 3:** Implementar. **Step 4:** PASS. Si la prueba 2 falla, seguir el spec §11 (no cambiar el criterio sin mostrar a Nick).
- [ ] **Step 5:** Commit `feat(motor): simulador, resultados y verificación contra legacy`.

### Task 9: dt adaptativo, convergencia y velocidad

**Files:** Modify `motor/simulador.py`; Create `tests/test_convergencia.py`

**Interfaces:** Sin cambios en la API. Reglas del spec §5 punto 5 con los valores del Global Constraints; en `control.modo_ssr = "onoff"`, dt ≤ t_ciclo.

- [ ] **Step 1: Tests:** `test_dt_mitad` (prueba 4) < 0.2 pp; `test_malla_doble` (prueba 5) < 0.5 pp y el test imprime el rango; `test_rampa_rapida_meseta_cero`: rampa 50, meseta 0 → termina en < 60 s; `test_velocidad` (prueba 6): caso base ≤ 5 s (marcar `lento`, documentar la máquina en el mensaje del assert).
- [ ] **Step 2–4:** FAIL → implementar → PASS. **Step 5:** Commit.

### Task 10: Lote, gráficas y README

**Files:** Create `motor/lote.py`, `motor/graficas.py`, `tests/test_lote.py`; Modify `motor/__init__.py` (exporta `simular_lote`), `README.md` (cómo correr un caso, un lote y las pruebas)

**Interfaces:** Produces `simular_lote(lista: list[Parametros], procesos: int = 12) -> list[Resultado | Exception]` (orden de entrada; `ProcessPoolExecutor`; variables `OMP/OPENBLAS/MKL/VECLIB_*_NUM_THREADS = "1"`); `graficar(resultado, carpeta) -> list[str]` (temperaturas vs. t, masa residual, mapa 2D final).

- [ ] **Step 1: Tests:** `test_lote_con_falla`: 3 casos, uno que lanza (monkeypatch) → lista de 3 con una excepción en su posición; `test_graficas`: crea los PNG en `tmp_path`.
- [ ] **Step 2–4:** FAIL → implementar → PASS. **Step 5:** Commit.

### Task 11: Verificación en MATLAB (la corre Nick)

**Files:** Create `matlab/verificacion/exportar_caso.py` (escribe `caso_conduccion.csv` con malla, k, ρcp, T₀, T_horno y T(t) de Python), `matlab/verificacion/verificar_conduccion.m` (ode15s sobre la misma malla sin reacción, lee el CSV y reporta error máx), `matlab/verificacion/README.md`

- [ ] **Step 1:** `exportar_caso.py` corre y crea el CSV (test `tests/test_exportar_matlab.py`).
- [ ] **Step 2:** Nick corre `verificar_conduccion` en MATLAB y pega el número. Criterio: error ≤ 1 % de ΔT₀ (prueba 7).
- [ ] **Step 3:** Commit `test: verificación térmica en MATLAB`.

### Cierre

- [ ] `python -m pytest tests -q` (incluye `-m lento`) → todo PASS; pegar el resumen en la bitácora con los tiempos.
- [ ] Copiar la rama al Mac (`modelos/PyrolysisReactorCoffee`) y correr el caso base en la Terminal de macOS.
- [ ] Entrada HECHO en la bitácora (GitHub + Drive) y anunciar que la interfaz quedó fija (desbloquea a Cursor/Copilot).
- [ ] PROPUESTA: borrar `legacy/` (solo con aprobación de Nick).
