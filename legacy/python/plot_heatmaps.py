"""
plot_heatmaps.py
Generador de gráficos de alta calidad para la visualización térmica y química:
1. Mapas de calor 2D del reactor cilíndrico (r vs z) en diferentes momentos del proceso
2. Perfiles dinámicos de temperatura (Horno, Pared, Centro con L'Hôpital)
3. Curva termogravimétrica de masa residual (TGA simulada)
4. Diagrama de barras de rendimientos de productos (Biochar, Biogás, Biovolátiles)
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_all_results(results, output_dir="."):
    """Genera y guarda todas las figuras analíticas del reactor"""
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    # --------------------------------------------------------------------------
    # 1. MAPA DE CALOR 2D EVOLUTIVO (r vs z)
    # --------------------------------------------------------------------------
    snapshots = results['snapshots']
    n_snaps = len(snapshots)
    geo = results['geo']
    
    # Malla de coordenadas en mm y cm
    R_grid_mm, Z_grid_cm = np.meshgrid(geo.r_grid * 1000.0, geo.z_grid * 100.0, indexing='ij')
    
    fig, axes = plt.subplots(1, n_snaps, figsize=(3.8 * n_snaps, 6), sharey=True)
    if n_snaps == 1:
        axes = [axes]
        
    vmin = min([np.min(s['T_map'] - 273.15) for s in snapshots])
    vmax = max([np.max(s['T_map'] - 273.15) for s in snapshots])
    
    for idx, (ax, snap) in enumerate(zip(axes, snapshots)):
        T_degC = snap['T_map'] - 273.15
        
        # Mapa de contorno continuo
        cs = ax.contourf(
            R_grid_mm, Z_grid_cm, T_degC,
            levels=np.linspace(vmin, vmax, 40),
            cmap='inferno'
        )
        # Líneas de contorno con etiquetas
        lines = ax.contour(
            R_grid_mm, Z_grid_cm, T_degC,
            levels=np.linspace(vmin, vmax, 8),
            colors='white', linewidths=0.7, alpha=0.6
        )
        ax.clabel(lines, inline=True, fontsize=8, fmt='%.0f°C')
        
        ax.set_title(f"{snap['label']}", fontsize=11, fontweight='bold')
        ax.set_xlabel("Radio r (mm)", fontsize=10)
        if idx == 0:
            ax.set_ylabel("Altura axial z (cm)", fontsize=10)
            
    cbar = fig.colorbar(cs, ax=axes, orientation='horizontal', fraction=0.04, pad=0.15)
    cbar.set_label("Temperatura (°C)", fontsize=11, fontweight='bold')
    plt.suptitle("Evolución del Mapa Térmico 2D en el Lecho Vertical de Residuos de Café", fontsize=13, fontweight='bold', y=0.98)
    
    map_filename = f"{output_dir}/mapa_de_calor_2d.png"
    plt.savefig(map_filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Guardado mapa de calor 2D: {map_filename}")

    # --------------------------------------------------------------------------
    # 2. PERFILES DINÁMICOS DE TEMPERATURA Y CURVAS DE RENDIMIENTO
    # --------------------------------------------------------------------------
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    t_min = results['time_min']
    
    # Panel 1: Curvas de Temperatura en el Tiempo
    ax1.plot(t_min, results['T_setpoint'], 'k--', label="Setpoint Programado (PID)", linewidth=1.5)
    ax1.plot(t_min, results['T_furnace'], 'm-', label="Horno (Inercia piedra refractaria)", linewidth=2.0)
    ax1.plot(t_min, results['T_wall'], 'r-', label="Pared exterior reactor (r = R)", linewidth=2.0)
    ax1.plot(t_min, results['T_mid'], 'g-.', label="Posición intermedia (r = R/2)", linewidth=1.8)
    ax1.plot(t_min, results['T_center'], 'b-', label="Centro del lecho (r = 0, L'Hôpital)", linewidth=2.2)
    ax1.set_xlabel("Tiempo (minutos)", fontsize=11)
    ax1.set_ylabel("Temperatura (°C)", fontsize=11)
    ax1.set_title("Evolución Térmica en el Lecho Vertical vs Horno", fontsize=12, fontweight='bold')
    ax1.legend(loc='lower right', frameon=True)
    ax1.grid(True, linestyle=':', alpha=0.6)
    
    # Panel 2: Perfiles Radiales T(r) en distintas alturas
    # Tomar la última instantánea
    last_snap = snapshots[-1]
    mid_z = geo.Nz // 2
    ax2.plot(geo.r_grid * 1000.0, last_snap['T_map'][:, 0] - 273.15, 'c-o', label="Entrada superior (z = 0)", markersize=4)
    ax2.plot(geo.r_grid * 1000.0, last_snap['T_map'][:, mid_z] - 273.15, 'b-s', label="Zona media (z = L/2)", markersize=4)
    ax2.plot(geo.r_grid * 1000.0, last_snap['T_map'][:, -1] - 273.15, 'r-^', label="Salida inferior (z = L)", markersize=4)
    ax2.set_xlabel("Radio r (mm)", fontsize=11)
    ax2.set_ylabel("Temperatura (°C)", fontsize=11)
    ax2.set_title("Perfil Radial de Temperatura T(r) al Final del Proceso", fontsize=12, fontweight='bold')
    ax2.legend(loc='best', frameon=True)
    ax2.grid(True, linestyle=':', alpha=0.6)
    
    # Panel 3: Curva Termogravimétrica (Pérdida de masa residual)
    ax3.plot(t_min, results['residual_mass_pct'], 'k-', linewidth=2.5, label="Masa sólida residual (% wt)")
    ax3.set_xlabel("Tiempo (minutos)", fontsize=11)
    ax3.set_ylabel("Masa residual (% wt)", fontsize=11)
    ax3.set_title("Curva de Pérdida de Masa del Café (Termogravimetría en Horno)", fontsize=12, fontweight='bold')
    ax3.set_ylim([0, 105])
    ax3.legend(loc='upper right', frameon=True)
    ax3.grid(True, linestyle=':', alpha=0.6)
    
    # Panel 4: Evolución Acumulada de Rendimientos
    ax4.plot(t_min, results['yield_biochar'], color='#4a3525', linewidth=2.2, label="Biochar")
    ax4.plot(t_min, results['yield_biovolatiles'], color='#d95f02', linewidth=2.2, label="Biovolátiles / Bio-oil")
    ax4.plot(t_min, results['yield_biogas'], color='#7570b3', linewidth=2.2, label="Biogás")
    ax4.plot(t_min, results['yield_water'], color='#1b9e77', linewidth=1.8, linestyle='--', label="Vapor H2O")
    ax4.set_xlabel("Tiempo (minutos)", fontsize=11)
    ax4.set_ylabel("Rendimiento acumulado (% wt)", fontsize=11)
    ax4.set_title("Evolución de Rendimientos de Productos de Pirólisis", fontsize=12, fontweight='bold')
    ax4.legend(loc='center left', frameon=True)
    ax4.grid(True, linestyle=':', alpha=0.6)
    
    plt.tight_layout()
    profiles_filename = f"{output_dir}/perfiles_temperatura_y_masa.png"
    plt.savefig(profiles_filename, dpi=300)
    plt.close()
    print(f"  [OK] Guardado panel de perfiles y masa: {profiles_filename}")

    # --------------------------------------------------------------------------
    # 3. DIAGRAMA DE BARRAS DE RENDIMIENTOS FINALES
    # --------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5.5))
    yields = results['final_yields']
    
    categories = ['Biochar\n(Sólido)', 'Biovolátiles / Bio-oil\n(Alquitranes + Condensables)', 'Biogás\n(No condensables)', 'Vapor H2O\n(Humedad + Reacción)']
    values = [
        yields['biochar_pct'],
        yields['biovolatiles_pct'],
        yields['biogas_pct'],
        yields['water_vapor_pct']
    ]
    colors = ['#3E2723', '#BF360C', '#1A237E', '#00695C']
    
    bars = ax.bar(categories, values, color=colors, width=0.55, edgecolor='black', linewidth=1.2)
    ax.set_ylabel("Rendimiento (% en peso)", fontsize=12, fontweight='bold')
    ax.set_title(
        f"Distribución de Productos de Pirólisis de Residuos de Café\n(Rampa {results['ramp'].ramp_rate_K_s*60:.0f}°C/min a {(results['ramp'].T_final_K-273.15):.0f}°C)",
        fontsize=12, fontweight='bold'
    )
    ax.set_ylim([0, max(values) * 1.25])
    
    for bar in bars:
        height = bar.get_height()
        ax.annotate(
            f"{height:.1f}%",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 4),
            textcoords="offset points",
            ha='center', va='bottom',
            fontsize=11, fontweight='bold'
        )
        
    ax.grid(axis='y', linestyle=':', alpha=0.7)
    plt.tight_layout()
    bar_filename = f"{output_dir}/rendimientos_finales_barras.png"
    plt.savefig(bar_filename, dpi=300)
    plt.close()
    print(f"  [OK] Guardado gráfico de rendimientos finales: {bar_filename}")
