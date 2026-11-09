"""
plot_lumped.py
Generador de figuras analíticas para el modelo cinético agrupado (Lumped Model):
1. Mapa de calor 2D evolutivo (r vs z)
2. Perfiles térmicos y panel explícito de CINÉTICAS DE CONSUMO Y PRODUCCIÓN:
   - Tasa de consumo de biomasa de café (r_biomass)
   - Tasa de formación de biovolátiles (r_volatiles)
   - Tasa de formación de biochar (r_char)
   - Tasa de formación de Gases No Condensables (r_gnc)
3. Curva de pérdida de masa residual
4. Diagrama de barras de rendimientos finales con Gases No Condensables (GNC)
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_lumped_results(results, output_dir="."):
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    geo = results['geo']
    t_min = results['time_min']
    snapshots = results['snapshots']
    
    # --------------------------------------------------------------------------
    # 1. MAPA DE CALOR 2D
    # --------------------------------------------------------------------------
    n_snaps = len(snapshots)
    R_grid_mm, Z_grid_cm = np.meshgrid(geo.r_grid * 1000.0, geo.z_grid * 100.0, indexing='ij')
    
    fig, axes = plt.subplots(1, n_snaps, figsize=(3.8 * n_snaps, 6), sharey=True)
    if n_snaps == 1:
        axes = [axes]
        
    vmin = min([np.min(s['T_map'] - 273.15) for s in snapshots])
    vmax = max([np.max(s['T_map'] - 273.15) for s in snapshots])
    
    for idx, (ax, snap) in enumerate(zip(axes, snapshots)):
        T_degC = snap['T_map'] - 273.15
        cs = ax.contourf(R_grid_mm, Z_grid_cm, T_degC, levels=np.linspace(vmin, vmax, 40), cmap='inferno')
        lines = ax.contour(R_grid_mm, Z_grid_cm, T_degC, levels=np.linspace(vmin, vmax, 8), colors='white', linewidths=0.7, alpha=0.6)
        ax.clabel(lines, inline=True, fontsize=8, fmt='%.0f°C')
        ax.set_title(f"{snap['label']}", fontsize=11, fontweight='bold')
        ax.set_xlabel("Radio r (mm)", fontsize=10)
        if idx == 0:
            ax.set_ylabel("Altura axial z (cm)", fontsize=10)
            
    cbar = fig.colorbar(cs, ax=axes, orientation='horizontal', fraction=0.04, pad=0.15)
    cbar.set_label("Temperatura (°C)", fontsize=11, fontweight='bold')
    plt.suptitle("Mapa Térmico 2D en Horno Vertical (Modelo Agrupado - Lumped)", fontsize=13, fontweight='bold', y=0.98)
    
    map_file = f"{output_dir}/mapa_de_calor_lumped.png"
    plt.savefig(map_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Guardado mapa de calor 2D (Lumped): {map_file}")

    # --------------------------------------------------------------------------
    # 2. PERFILES TÉRMICOS, CINÉTICAS DE CONSUMO Y CURVAS DE MASA
    # --------------------------------------------------------------------------
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10.5))
    
    # Panel 1: Temperaturas en el tiempo
    ax1.plot(t_min, results['T_setpoint'], 'k--', label="Setpoint Programado (PID)", linewidth=1.5)
    ax1.plot(t_min, results['T_furnace'], 'm-', label="Horno (Inercia piedra refractaria)", linewidth=2.0)
    ax1.plot(t_min, results['T_wall'], 'r-', label="Pared exterior (r = R)", linewidth=2.0)
    ax1.plot(t_min, results['T_mid'], 'g-.', label="Posición media (r = R/2)", linewidth=1.8)
    ax1.plot(t_min, results['T_center'], 'b-', label="Centro (r = 0, L'Hôpital)", linewidth=2.2)
    ax1.set_xlabel("Tiempo (minutos)", fontsize=11)
    ax1.set_ylabel("Temperatura (°C)", fontsize=11)
    ax1.set_title("Evolución de Temperaturas en el Lecho", fontsize=12, fontweight='bold')
    ax1.legend(loc='lower right', frameon=True)
    ax1.grid(True, linestyle=':', alpha=0.6)
    
    # Panel 2: TASAS CINÉTICAS EXPLÍCITAS (Consumo de Biomasa vs Formación de Volátiles y Char)
    ax2.plot(t_min, results['rate_biomass_decay'], 'k-', linewidth=2.2, label="Consumo Biomasa Café (-dρB/dt)")
    ax2.plot(t_min, results['rate_volatiles'], color='#d95f02', linewidth=2.0, label="Formación Biovolátiles (r_vol)")
    ax2.plot(t_min, results['rate_char'], color='#4a3525', linewidth=2.0, label="Formación Biochar (r_char)")
    ax2.plot(t_min, results['rate_gnc'], color='#1A237E', linewidth=1.8, linestyle='--', label="Formación Gases No Condensables (r_gnc)")
    ax2.set_xlabel("Tiempo (minutos)", fontsize=11)
    ax2.set_ylabel("Velocidad de Reacción (kg / (m³·h))", fontsize=11)
    ax2.set_title("Cinética de Consumo de Biomasa y Producción de Fases", fontsize=12, fontweight='bold')
    ax2.legend(loc='upper right', frameon=True)
    ax2.grid(True, linestyle=':', alpha=0.6)
    
    # Panel 3: Curva Termogravimétrica (Pérdida de masa)
    ax3.plot(t_min, results['residual_mass_pct'], 'k-', linewidth=2.5, label="Masa sólida residual (% wt)")
    ax3.set_xlabel("Tiempo (minutos)", fontsize=11)
    ax3.set_ylabel("Masa sólida residual (% peso)", fontsize=11)
    ax3.set_title("Curva de Pérdida de Masa (TGA Simulada)", fontsize=12, fontweight='bold')
    ax3.set_ylim([0, 105])
    ax3.legend(loc='upper right', frameon=True)
    ax3.grid(True, linestyle=':', alpha=0.6)
    
    # Panel 4: Rendimientos acumulados en el tiempo
    ax4.plot(t_min, results['yield_biochar'], color='#4a3525', linewidth=2.2, label="Biochar")
    ax4.plot(t_min, results['yield_biovolatiles'], color='#d95f02', linewidth=2.2, label="Biovolátiles / Bio-oil")
    ax4.plot(t_min, results['yield_gnc'], color='#1A237E', linewidth=2.2, label="Gases No Condensables (GNC)")
    ax4.plot(t_min, results['yield_water'], color='#1b9e77', linewidth=1.8, linestyle='--', label="Vapor H2O")
    ax4.set_xlabel("Tiempo (minutos)", fontsize=11)
    ax4.set_ylabel("Rendimiento acumulado (% wt)", fontsize=11)
    ax4.set_title("Evolución de Rendimientos Acumulados (Lumped)", fontsize=12, fontweight='bold')
    ax4.legend(loc='center left', frameon=True)
    ax4.grid(True, linestyle=':', alpha=0.6)
    
    plt.tight_layout()
    prof_file = f"{output_dir}/perfiles_y_cineticas_lumped.png"
    plt.savefig(prof_file, dpi=300)
    plt.close()
    print(f"  [OK] Guardado panel de perfiles y cinéticas (Lumped): {prof_file}")

    # --------------------------------------------------------------------------
    # 3. DIAGRAMA DE BARRAS DE RENDIMIENTOS FINALES
    # --------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    yields = results['final_yields']
    
    categories = [
        'Biochar\n(Sólido)',
        'Biovolátiles / Bio-oil\n(Condensables)',
        'Gases No Condensables\n(GNC: CO, CO2, CH4, H2)',
        'Vapor H2O\n(Humedad inicial)'
    ]
    values = [
        yields['biochar_pct'],
        yields['biovolatiles_pct'],
        yields['gnc_pct'],
        yields['water_vapor_pct']
    ]
    colors = ['#3E2723', '#BF360C', '#1A237E', '#00695C']
    
    bars = ax.bar(categories, values, color=colors, width=0.52, edgecolor='black', linewidth=1.2)
    ax.set_ylabel("Rendimiento (% en peso)", fontsize=12, fontweight='bold')
    ax.set_title(
        f"Distribución Final de Productos - Modelo Cinético Agrupado (Lumped)\n(Rampa {results['ramp'].ramp_rate_K_s*60:.0f}°C/min hasta {(results['ramp'].T_final_K-273.15):.0f}°C)",
        fontsize=12, fontweight='bold'
    )
    ax.set_ylim([0, max(values) * 1.25])
    
    for bar in bars:
        h = bar.get_height()
        ax.annotate(
            f"{h:.1f}%",
            xy=(bar.get_x() + bar.get_width() / 2, h),
            xytext=(0, 4),
            textcoords="offset points",
            ha='center', va='bottom',
            fontsize=11, fontweight='bold'
        )
        
    ax.grid(axis='y', linestyle=':', alpha=0.7)
    plt.tight_layout()
    bar_file = f"{output_dir}/rendimientos_lumped_barras.png"
    plt.savefig(bar_file, dpi=300)
    plt.close()
    print(f"  [OK] Guardado gráfico de barras (Lumped): {bar_file}")
