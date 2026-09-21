"""
run_lumped.py (Raíz del proyecto)
Launcher de conveniencia para el modelo cinético agrupado (LUMPED).
Permite ejecutar directamente desde la carpeta principal:
  python run_lumped.py --ramp 10 --temp 500 --soak 30
"""

import os
import sys
import subprocess

root_dir = os.path.dirname(os.path.abspath(__file__))
venv_python = os.path.join(root_dir, ".venv", "Scripts", "python.exe")

if os.path.exists(venv_python) and os.path.normcase(sys.executable) != os.path.normcase(venv_python):
    result = subprocess.run([venv_python] + sys.argv)
    sys.exit(result.returncode)

python_dir = os.path.join(root_dir, "python")
sys.path.insert(0, python_dir)
os.chdir(python_dir)

import run_lumped

if __name__ == "__main__":
    parser = run_lumped.argparse.ArgumentParser(description="Simulador de Pirólisis con Modelo Cinético Agrupado (Lumped)")
    parser.add_argument("--ramp", type=float, default=10.0, help="Tasa de calentamiento en °C/min (def: 10.0)")
    parser.add_argument("--temp", type=float, default=500.0, help="Temperatura máxima en °C (def: 500.0)")
    parser.add_argument("--soak", type=float, default=45.0, help="Tiempo de meseta en min (def: 45.0)")
    parser.add_argument("--outdir", type=str, default=".", help="Directorio de salida para figuras")
    
    args = parser.parse_args()
    run_lumped.run_lumped_simulation(ramp_rate=args.ramp, final_temp=args.temp, soak_time=args.soak, output_dir=args.outdir)
