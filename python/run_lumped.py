"""
run_lumped.py
Script ejecutable para correr la simulación con el modelo cinético agrupado (Lumped).

Uso:
  python run_lumped.py
  python run_lumped.py --ramp 10 --temp 500 --soak 30
  python run_lumped.py --ramp 5 --temp 450 --soak 45
"""

import sys
import argparse
from config import GeometryConfig, CoffeeBiomassConfig, ThermalPropertiesConfig, FurnaceConfig, HeatingRampConfig
from simulation_lumped import LumpedPyrolysisSimulator
from plot_lumped import plot_lumped_results

def run_lumped_simulation(ramp_rate=10.0, final_temp=500.0, soak_time=45.0, output_dir="."):
    print(f"\n>>> INICIANDO SIMULACIÓN LUMPED: RAMPA = {ramp_rate} °C/min, T_FINAL = {final_temp} °C, SOAK = {soak_time} min <<<")
    
    geo_cfg = GeometryConfig(R_inner=0.025, L_bed=0.20, Nr=15, Nz=20)
    bio_cfg = CoffeeBiomassConfig(moisture_fraction=0.08)
    therm_cfg = ThermalPropertiesConfig()
    furnace_cfg = FurnaceConfig(power_kW=2.5, stone_thermal_inertia=240.0)
    ramp_cfg = HeatingRampConfig(
        T_initial_C=25.0,
        ramp_rate_C_per_min=ramp_rate,
        T_final_C=final_temp,
        soak_time_min=soak_time
    )
    
    sim = LumpedPyrolysisSimulator(geo_cfg, bio_cfg, therm_cfg, furnace_cfg, ramp_cfg)
    results = sim.run(max_dt=0.5, save_interval_s=20.0, verbose=True)
    
    print("\nGenerando figuras y mapas de calor (Lumped)...")
    plot_lumped_results(results, output_dir=output_dir)
    print("\nSimulación Lumped completada con éxito.")
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulador de Pirólisis con Modelo Cinético Agrupado (Lumped)")
    parser.add_argument("--ramp", type=float, default=10.0, help="Tasa de calentamiento en °C/min (def: 10.0)")
    parser.add_argument("--temp", type=float, default=500.0, help="Temperatura máxima en °C (def: 500.0)")
    parser.add_argument("--soak", type=float, default=45.0, help="Tiempo de meseta en min (def: 45.0)")
    parser.add_argument("--outdir", type=str, default=".", help="Directorio de salida para figuras")
    
    args = parser.parse_args()
    run_lumped_simulation(ramp_rate=args.ramp, final_temp=args.temp, soak_time=args.soak, output_dir=args.outdir)
