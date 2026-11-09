# Spec — Ensayo térmico y calibración (subproyecto B)

Fecha: 2026-11-08 (Bogotá). Estado: **borrador para aprobar por Nick**. Base: spec del motor §6 y §11, `parametros/base.toml` (main, ed907cf).

## 1. Para qué es

Medir dos curvas reales (tubo vacío y tubo con café) y ajustar con ellas los parámetros térmicos del motor que hoy son [Supuesto] o [Por validar]. Con el modelo calibrado se sintoniza el PID (subproyecto C).

## 2. Qué se calibra y con qué corrida

| Parámetro (base.toml) | Valor hoy | Se calibra con | Señal que lo fija |
|---|---|---|---|
| `horno.C_h` | 20000 J/K [Por validar] | Corrida 1 (vacío) | Pendiente inicial de T_horno con potencia conocida: dT/dt = P/C_h |
| `horno.UA_p` | 5 W/K [Por validar] | Corrida 1 | Estado estable: P = UA_p·(T_horno − T_amb); enfriamiento: τ = C_h/UA_p |
| `horno.h_c` | 90 W/m²K [Supuesto] | Corrida 1 | Retraso T_horno − T_pared en la rampa |
| `horno.emisividad` | 0 (radiación dentro de h_c) | Corrida 1 | Cómo cambia ese retraso con la temperatura |
| `lecho.k_virgen_0` | 0.12 W/mK [Supuesto] | Corrida 2, tramo 25–200 °C | Retraso T_pared − T_centro antes de la pirólisis |
| `lecho.k_char_0` | 0.08 W/mK [Supuesto] | Corrida 2, enfriamiento final | Retraso T_pared − T_centro con el lecho ya carbonizado y sin reacciones |
| `lecho.R_contacto_extra` | 0.011 m²K/W [Calculado] | Corrida 2 (ajuste conjunto con k) | Débil con 2 termopares (ver §4) |
| `termopares.tau` | 10 s [Supuesto] | Prueba aparte (§3.1) | Escalón en agua hirviendo |
| `horno.P_max` | 2500 W [Supuesto] | Medición aparte (§3.1) | Pinza amperimétrica + voltaje |

No se calibran (quedan como están): pendientes `k_*_pend` y `cp_*` (no se separan con estas dos corridas), cinética y geometría.
Kp, Ki, Kd, t_ciclo **no** se calibran aquí: se sintonizan en C con el modelo ya calibrado.

## 3. Montaje y mediciones

### 3.1 Antes de las corridas (barato, una vez)
- **τ de cada termopar:** pasar de aire a agua hirviendo y medir el tiempo al 63 %. Reemplaza `termopares.tau` = 10 s.
- **P_max real:** corriente con pinza y voltaje de la red con la resistencia caliente y el SSR al 100 %. P = V·I.
- **Geometría real del tubo:** diámetro interno, espesor y largo calentado (hoy 25 mm / 2.5 mm / 0.20 m [Supuesto]).
- **Carga:** masa de café en balanza y altura del lecho → densidad aparente real. Humedad de una submuestra en estufa (hoy 8 % bh [Supuesto]). Referencia: el tubo base cabe ≈ 0.149 kg a 380 kg/m³ [Calculado: π·0.025²·0.20·380].

### 3.2 Termopares (tipo K con vaina, todos a z = 0.10 m, mitad del lecho)
1. **T_pared:** cara externa del tubo, pegado con cemento cerámico o abrazadera y tapado del lado del horno (si no, ve la radiación de la resistencia y mide de más). Es el mismo nodo que `historia.T_pared`.
2. **T_centro:** r = 0, entrando desde arriba por una vaina delgada. Mismo nodo que `historia.T_centro`.
3. **T_horno (muy recomendado):** aire de la cámara, entre resistencia y tubo, sin tocar ninguno. Sin él, C_h y UA_p solo se pueden sacar si el registro guarda la salida del SSR.
4. **Opcional, para separar R_contacto de k:** uno en el lecho cerca de la pared (r ≈ 20 mm).

### 3.3 Registro
Un CSV por corrida en `01_Horno_Vertical/datos/ensayo_termico/` (no va al repo), muestreo 1 s, columnas:
`t_s, T_horno_C, T_pared_C, T_centro_C, T_control_C, T_set_C, u_ssr_pct, T_amb_C`.
`u_ssr_pct` (salida del MC-5438) es clave: sin ella no se conoce la potencia. Depende del manual del MC-5438 (pendiente de Perplexity).
El mismo flujo de N₂ en las dos corridas (el gas también se lleva calor).

## 4. Las dos corridas

### Corrida 1 — tubo vacío (sin café, tapado igual que con carga)
Programa (si el MC-5438 permite salida manual):
1. SSR fijo a un % conocido (ej. 50 %) desde ambiente hasta estado estable (T casi constante 10 min).
2. Subir a 100 % hasta ≈ 550 °C.
3. Apagar y dejar enfriar registrando hasta < 100 °C.

Si no hay salida manual: rampa 10 °C/min a 550 °C con el PID, 20 min de meseta y enfriamiento libre, registrando `u_ssr_pct`.

**Qué da:** C_h y UA_p (calentamiento, estable y enfriamiento), y h_c y ε. Con el tubo vacío no hay lecho, así que k y R_contacto no aparecen: h_c queda fijado sin mezclarse con ellos.
Orden de magnitud con los valores de hoy [Calculado]: tubo de acero ≈ 0.65 kg → 326 J/K; h_c·A = 90 × 0.0346 m² = 3.1 W/K; τ del tubo ≈ 105 s, o sea que a 10 °C/min T_pared va ≈ 17 K detrás de T_horno. Ese retraso es lo que se mide.
**h_c contra ε:** con ε ≈ 0.7 [Por validar, Incropera A.11] la parte radiativa vale ≈ 8 W/m²K a 100 °C y ≈ 73 W/m²K a 500 °C [Calculado]. Si el retraso medido baja con la temperatura, la radiación pesa y se ajustan h_c y ε por separado; si no cambia, se deja ε = 0 y un solo h_c.

### Corrida 2 — con café
Caso base: 10 °C/min hasta 500 °C, 30 min de meseta, luego enfriamiento libre hasta < 100 °C registrando todo. C_h, UA_p, h_c y ε quedan **fijos** de la corrida 1.
Tramos que se usan:
- **25–200 °C (antes de pirolizar):** k_virgen_0 (+ R_contacto). Ojo con la meseta de secado cerca de 100 °C: por eso hay que medir la humedad real.
- **Pirólisis (200–500 °C):** no se ajusta; solo se compara (aquí pesan la cinética y el calor de reacción).
- **Enfriamiento final:** lecho ya carbonizado y sin reacciones → k_char_0 limpio.
Al final pesar el biochar: comparación gruesa con el 24.22 % del caso base.

### Por qué h_c y k del lecho no se confunden (riesgo §11)
- h_c sale solo de la corrida 1, donde no hay lecho. En la 2 queda fijo.
- En la corrida 2 el número de Biot del lecho es ≈ 9 [Calculado: h_global ≈ 45 W/m²K × 0.025 m / 0.12 W/mK]: el retraso del centro lo manda la conducción dentro del lecho, así que T_centro sí fija k.
- R_contacto (0.011 m²K/W) es solo ≈ 10 % de la resistencia del lecho (R/2k ≈ 0.10 m²K/W) [Calculado]: con dos termopares queda débil y correlacionada con k. Se reporta la correlación del ajuste; si es alta, se deja R_contacto en 0.011 o se pone el termopar opcional (§3.2.4).

## 5. Cómo se ajusta
- Script en `calibracion/` (rama `ensayo-termico`), solo con la interfaz pública: `cargar_parametros`, `simular`, `simular_lote`.
- Mínimos cuadrados (scipy `least_squares`) sobre RMS de T_pared y T_centro en los tramos de §4. Corridas pesadas en la Terminal de macOS con el `.venv`.
- Corrida 1: el motor necesita un lecho. Dos opciones (Nick decide, ver PROPUESTA en la bitácora):
  - (a) Modelo 0D del tubo en el script de calibración (C_tubo·dT_pared/dt = h_ext·A·(T_horno − T_pared), con el mismo h_ext del motor). No toca el motor.
  - (b) Agregar al motor un modo "tubo vacío".
- Corrida 2: mejor con T_horno medido como entrada (separa horno de tubo). Hoy el motor solo tiene `horno.modo` = ideal o físico. Opción: agregar `horno.modo = "medido"` (lee T_horno(t) del CSV). También va como PROPUESTA.

## 6. Criterio de cierre de B (propuesta, Nick aprueba)
- Corrida 1: RMS T_pared ≤ 5 K y RMS T_horno ≤ 5 K.
- Corrida 2: RMS T_pared ≤ 5 K; RMS T_centro ≤ 10 K fuera del tramo de pirólisis y ≤ 20 K en él.
- Todos los valores calibrados pasan a `base.toml` con etiqueta [Medido/Calibrado], fuente = archivo de la corrida y rango del ajuste.
- Validación (opcional, recomendada): una corrida 3 a 20 °C/min que no se usa para ajustar. Responde además cuánto pesa la rampa en el horno real.

## 7. Depende de
- Manual del MC-5438: salida manual del SSR, registro de la salida, t_ciclo (pendiente Perplexity → Nick).
- Ubicación del termopar de control real (pregunta abierta).
- Contar con 3 termopares tipo K y un registrador de 1 s.
