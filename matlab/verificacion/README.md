# Verificación térmica en MATLAB (prueba 7 del spec)

Compara la conducción implícita del motor (Python) con el mismo problema resuelto en MATLAB con `ode15s`.
El script de MATLAB vuelve a armar los volúmenes y las conductancias desde la geometría: no copia la matriz de Python.

Caso: malla base (14 + 4 nodos radiales, 20 axiales), lecho virgen con k y ρcp constantes, horno fijo a 500 °C,
h = 90 W/m²K, 30 min, sin reacción. Criterio: error máximo ≤ 1 % de (T_horno − T0) = 4.75 K.

1. En la Terminal, desde la raíz del repo: `.venv/bin/python matlab/verificacion/exportar_caso.py`
   (escribe `caso_conduccion.csv` y `caso_conduccion_meta.csv` aquí; también imprime el segundo camino en Python con BDF).
2. En MATLAB: `cd matlab/verificacion` y `verificar_conduccion`.
3. Pegar en la bitácora las dos últimas líneas que imprime.

Resultado de referencia (GNU Octave en Linux, 2026-11-08): error máximo 0.349 K → PASA. El BDF de Python da lo mismo (0.349 K).
