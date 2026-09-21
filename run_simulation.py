"""
run_simulation.py (Raíz del proyecto)
Launcher de conveniencia. Permite ejecutar la simulación directamente desde la carpeta raíz:
  python run_simulation.py --ramp 10 --temp 500 --soak 30

Detecta automáticamente el entorno virtual .venv y delega a python/run_simulation.py.
"""

import os
import sys
import subprocess

# 1. Si se ejecutó con el python global del sistema, re-ejecutar con el python del .venv
root_dir = os.path.dirname(os.path.abspath(__file__))
venv_python = os.path.join(root_dir, ".venv", "Scripts", "python.exe")

if os.path.exists(venv_python) and os.path.normcase(sys.executable) != os.path.normcase(venv_python):
    result = subprocess.run([venv_python] + sys.argv)
    sys.exit(result.returncode)

# 2. Agregar carpeta python/ al path y cambiar directorio de trabajo a python/
python_dir = os.path.join(root_dir, "python")
sys.path.insert(0, python_dir)
os.chdir(python_dir)

# 3. Importar y ejecutar el script principal
import run_simulation

if __name__ == "__main__":
    parser = run_simulation.argparse.ArgumentParser(description="Simulador de Pirólisis de Residuos de Café en Horno Vertical")
    parser.add_argument("--ramp", type=float, default=10.0, help="Tasa de calentamiento en °C/min (def: 10.0)")
    parser.add_argument("--temp", type=float, default=500.0, help="Temperatura máxima en °C (def: 500.0)")
    parser.add_argument("--soak", type=float, default=45.0, help="Tiempo de meseta isotérmica en min (def: 45.0)")
    parser.add_argument("--compare", action="store_true", help="Ejecutar comparativa entre rampa lenta (5°C/min) y rápida (20°C/min)")
    parser.add_argument("--outdir", type=str, default=".", help="Directorio de salida para gráficos")
    
    args = parser.parse_args()
    if args.compare:
        run_simulation.compare_ramps(output_dir=args.outdir)
    else:
        run_simulation.run_single_simulation(ramp_rate=args.ramp, final_temp=args.temp, soak_time=args.soak, output_dir=args.outdir)
