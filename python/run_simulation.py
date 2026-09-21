"""
run_simulation.py
Script ejecutable principal para el reactor de pirólisis de residuos de café.
Permite ejecutar simulaciones individuales o comparar diferentes rampas de temperatura.

Uso:
  python run_simulation.py
  python run_simulation.py --ramp 5 --temp 450 --soak 30
  python run_simulation.py --compare
"""

import sys
import argparse
from config import GeometryConfig, CoffeeBiomassConfig, ThermalPropertiesConfig, FurnaceConfig, HeatingRampConfig
from simulation import PyrolysisSimulator
from plot_heatmaps import plot_all_results
import matplotlib.pyplot as plt
import numpy as np

def run_single_simulation(ramp_rate=10.0, final_temp=500.0, soak_time=45.0, output_dir="."):
    print(f"\n>>> INICIANDO SIMULACIÓN CON RAMPA = {ramp_rate} °C/min, T_FINAL = {final_temp} °C, SOAK = {soak_time} min <<<")
    
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
    
    sim = PyrolysisSimulator(geo_cfg, bio_cfg, therm_cfg, furnace_cfg, ramp_cfg)
    results = sim.run(max_dt=0.5, save_interval_s=20.0, verbose=True)
    
    print("\nGenerando figuras y mapas de calor...")
    plot_all_results(results, output_dir=output_dir)
    print("\nSimulación completada con éxito.")
    return results

def compare_ramps(output_dir="."):
    """Compara una rampa lenta (5 °C/min) vs una rampa rápida (20 °C/min)"""
    print("\n==================================================================")
    print("  COMPARACIÓN DE RAMPAS TÉRMICAS: 5 °C/min vs 20 °C/min")
    print("==================================================================")
    
    # Caso 1: Rampa lenta (Favorece Biochar)
    ramp_slow = HeatingRampConfig(T_initial_C=25.0, ramp_rate_C_per_min=5.0, T_final_C=500.0, soak_time_min=30.0)
    sim1 = PyrolysisSimulator(ramp_cfg=ramp_slow)
    res1 = sim1.run(max_dt=0.5, save_interval_s=30.0, verbose=False)
    
    # Caso 2: Rampa rápida (Favorece Biovolátiles / Bio-oil)
    ramp_fast = HeatingRampConfig(T_initial_C=25.0, ramp_rate_C_per_min=20.0, T_final_C=500.0, soak_time_min=30.0)
    sim2 = PyrolysisSimulator(ramp_cfg=ramp_fast)
    res2 = sim2.run(max_dt=0.5, save_interval_s=30.0, verbose=False)
    
    # Graficar comparación
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    
    # Panel 1: Curvas de masa residual
    ax1.plot(res1['time_min'], res1['residual_mass_pct'], 'b-', linewidth=2.2, label="Rampa Lenta (5 °C/min)")
    ax1.plot(res2['time_min'], res2['residual_mass_pct'], 'r--', linewidth=2.2, label="Rampa Rápida (20 °C/min)")
    ax1.set_xlabel("Tiempo (minutos)", fontsize=11)
    ax1.set_ylabel("Masa sólida residual (% wt)", fontsize=11)
    ax1.set_title("Efecto de la Rampa en la Velocidad de Pérdida de Masa", fontsize=12, fontweight='bold')
    ax1.legend(frameon=True)
    ax1.grid(True, linestyle=':', alpha=0.6)
    
    # Panel 2: Comparación de rendimientos finales
    categories = ['Biochar', 'Biovolátiles\n(Bio-oil)', 'Biogás', 'Vapor H2O']
    vals1 = [res1['final_yields']['biochar_pct'], res1['final_yields']['biovolatiles_pct'], res1['final_yields']['biogas_pct'], res1['final_yields']['water_vapor_pct']]
    vals2 = [res2['final_yields']['biochar_pct'], res2['final_yields']['biovolatiles_pct'], res2['final_yields']['biogas_pct'], res2['final_yields']['water_vapor_pct']]
    
    x = np.arange(len(categories))
    width = 0.35
    
    rects1 = ax2.bar(x - width/2, vals1, width, label='5 °C/min (Lenta)', color='#2b5c8f', edgecolor='black')
    rects2 = ax2.bar(x + width/2, vals2, width, label='20 °C/min (Rápida)', color='#d95f02', edgecolor='black')
    
    ax2.set_ylabel("Rendimiento (% peso)", fontsize=11, fontweight='bold')
    ax2.set_title("Distribución de Productos según la Rampa", fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(categories, fontsize=10)
    ax2.legend(frameon=True)
    ax2.grid(axis='y', linestyle=':', alpha=0.6)
    
    # Etiquetas sobre las barras
    for rect in rects1:
        h = rect.get_height()
        ax2.annotate(f"{h:.1f}%", xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9, fontweight='bold')
    for rect in rects2:
        h = rect.get_height()
        ax2.annotate(f"{h:.1f}%", xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9, fontweight='bold')
        
    plt.tight_layout()
    cmp_filename = f"{output_dir}/comparacion_de_rampas.png"
    plt.savefig(cmp_filename, dpi=300)
    plt.close()
    print(f"  [OK] Guardado gráfico comparativo de rampas: {cmp_filename}")
    
    print("\nResumen comparativo:")
    print(f"  • Rampa 5 °C/min:  Biochar = {vals1[0]:.1f}%, Bio-oil = {vals1[1]:.1f}%, Biogás = {vals1[2]:.1f}%")
    print(f"  • Rampa 20 °C/min: Biochar = {vals2[0]:.1f}%, Bio-oil = {vals2[1]:.1f}%, Biogás = {vals2[2]:.1f}%")
    print("  -> Conclusión: Rampas lentas favorecen la carbonización secundaria aumentando el rendimiento de Biochar.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulador de Pirólisis de Residuos de Café en Horno Vertical")
    parser.add_argument("--ramp", type=float, default=10.0, help="Tasa de calentamiento en °C/min (def: 10.0)")
    parser.add_argument("--temp", type=float, default=500.0, help="Temperatura máxima en °C (def: 500.0)")
    parser.add_argument("--soak", type=float, default=45.0, help="Tiempo de meseta isotérmica en min (def: 45.0)")
    parser.add_argument("--compare", action="store_true", help="Ejecutar comparativa entre rampa lenta (5°C/min) y rápida (20°C/min)")
    parser.add_argument("--outdir", type=str, default=".", help="Directorio de salida para gráficos")
    
    args = parser.parse_args()
    if args.compare:
        compare_ramps(output_dir=args.outdir)
    else:
        run_single_simulation(ramp_rate=args.ramp, final_temp=args.temp, soak_time=args.soak, output_dir=args.outdir)
