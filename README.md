# Simulador de Transferencia de Calor y Pirólisis de Café en Horno Vertical

Este proyecto contiene la implementación computacional completa en **Python** y en **MATLAB** del modelo matemático de transferencia de calor bidimensional axisimétrico $(r, z)$ acoplado a la cinética multicomponente de pirólisis de residuos de café (*Spent Coffee Grounds* - SCG).

El sistema modela:
1. **Horno eléctrico exterior**: Resistencias calefactoras, relé de estado sólido (SSR), controlador PID MaxThermo MC-5438 y recubrimiento refractario de piedra con alta inercia térmica.
2. **Reactor vertical**: Tubo de acero inoxidable con canastilla concéntrica de lecho poroso de biomasa, inyección superior de gas inerte ($N_2$) y descarga inferior de condensables/bio-oil ($\dot{m}_{oil}$).
3. **Ecuación de calor no lineal en coordenadas cilíndricas**:
   $$\rho_{eff}(T, X) \cdot C_{p,eff}(T, X) \frac{\partial T}{\partial t} = \frac{1}{r}\frac{\partial}{\partial r}\left( r \cdot k_{eff,r} \frac{\partial T}{\partial r} \right) + \frac{\partial}{\partial z}\left( k_{eff,z} \frac{\partial T}{\partial z} \right) - \rho_g C_{p,g} v_{g,z} \frac{\partial T}{\partial z} + \sum (-\Delta H_{r,i}) r_i$$
4. **Condición en el centro ($r=0$) mediante la Regla de L'Hôpital**:
   $$\lim_{r \to 0} \frac{1}{r}\frac{\partial T}{\partial r} = \left.\frac{\partial^2 T}{\partial r^2}\right|_{r=0} \implies \left.\frac{1}{r}\frac{\partial}{\partial r}\left(r k \frac{\partial T}{\partial r}\right)\right|_{r=0} = 2 k \left.\frac{\partial^2 T}{\partial r^2}\right|_{r=0} \approx 4 k \frac{T_1 - T_0}{\Delta r^2}$$
5. **Caracterización experimental de la biomasa**:
   Integración directa de los resultados de laboratorio del **GIEM (Universidad de Antioquia / EIA)** y espectroscopía **FTIR** (humedad 8% acondicionado, cenizas 2.31%, lípidos/aceite de café 15%, celulosa 12%, hemicelulosa 35%, lignina 28%).

---

## Estructura del Proyecto

```
pyrolysis_reactor/
├── python/
│   ├── config.py             # Parámetros geométricos, cinéticos, termofísicos y rampas
│   ├── kinetics.py           # Mecanismo cinético multicomponente adaptado a SCG
│   ├── furnace_control.py    # Horno con resistencias, piedra refractaria y PID MC-5438
│   ├── thermal_model.py      # Discretización espacial 2D (r, z) y L'Hôpital en r = 0
│   ├── simulation.py         # Motor de integración temporal y balances de masa
│   ├── plot_heatmaps.py      # Generador de mapas de calor 2D y curvas de rendimiento
│   └── run_simulation.py     # Script principal ejecutable con CLI
│
├── matlab/
│   ├── param.m               # Parámetros y propiedades físicas
│   ├── funciones_ct.m        # Sistema de EDOs para ode15s con L'Hôpital en r = 0
│   ├── simular_rampa.m       # Script principal de ejecución
│   └── graficar_resultados.m # Visualización 2D (contourf), perfiles y rendimientos
│
├── requirements.txt          # Dependencias de Python (numpy, scipy, matplotlib)
└── README.md
```

---

## ¿Cómo Ejecutar en Python?

El entorno virtual `.venv` ya se encuentra configurado en la carpeta del proyecto.

### 1. Ejecución Básica (Rampa por defecto: 10 °C/min hasta 500 °C, 45 min sostenimiento)
```powershell
cd C:\Users\Usuario\Documents\pyrolysis_reactor\python
..\.venv\Scripts\python.exe run_simulation.py
```

### 2. Modificar la Rampa y el Tiempo de Pirólisis
Puedes ingresar cualquier parámetro desde la línea de comandos:
```powershell
# Ejemplo: Rampa lenta de 5 °C/min hasta 450 °C, sostenida 30 min (mayor rendimiento de biochar)
..\.venv\Scripts\python.exe run_simulation.py --ramp 5 --temp 450 --soak 30

# Ejemplo: Rampa rápida de 20 °C/min hasta 550 °C, sostenida 20 min (mayor rendimiento de bio-oil)
..\.venv\Scripts\python.exe run_simulation.py --ramp 20 --temp 550 --soak 20
```

### 3. Comparativa Automática entre Rampas
```powershell
..\.venv\Scripts\python.exe run_simulation.py --compare
```

---

## ¿Cómo Ejecutar en MATLAB?

1. Abre MATLAB y navega a la carpeta:
   `C:\Users\Usuario\Documents\pyrolysis_reactor\matlab`
2. Para cambiar la rampa de calentamiento, abre `param.m` y modifica las líneas:
   ```matlab
   beta_C_min = 10.0;     % Tasa de calentamiento (°C/min)
   T_fin_C    = 500.0;    % Temperatura máxima isotérmica (°C)
   t_soak_min = 30.0;     % Tiempo de sostenimiento isotérmico (min)
   ```
3. Ejecuta en la consola de comandos de MATLAB:
   ```matlab
   simular_rampa
   ```
4. Se abrirán automáticamente 4 figuras interactivas:
   - **Figura 1**: Mapa de contorno 2D de calor (`contourf`) a lo largo del lecho ($r$ vs $z$).
   - **Figura 2**: Curvas de temperatura en el tiempo ($r=0$, $r=R/2$, pared $r=R$, horno).
   - **Figura 3**: Curva termogravimétrica de pérdida de masa ($W/W_0$).
   - **Figura 4**: Gráfico de barras con los porcentajes de **Biochar**, **Biovolátiles / Bio-oil**, **Biogás** y **Vapor de agua**.
