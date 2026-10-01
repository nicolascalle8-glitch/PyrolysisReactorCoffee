"""
simulador_interactivo.py
Simulador Completo e Interactivo de Transferencia de Calor y Pirólisis de Café.
CARACTERÍSTICAS:
1. Malla conjugada multicapa: Lecho de Café (0 a R_int) + Reactor de Acero Inoxidable (R_int a R_ext).
2. Cinética calibrada con sensibilidad drástica de rendimientos a la rampa (Lenta = Biochar, Rápida = Volátiles).
3. Panel científico de 6 gráficas (Especies vs Tiempo y vs Temperatura, TGA+DTG, Perfiles).
4. DOS MAPAS TÉRMICOS 2D INTERACTIVOS CON SLIDER:
   - Mapa Longitudinal 2D (r vs z) con delimitación explícita del tubo de acero inox.
   - Corte Transversal Circular 360° mostrando el anillo de acero inox envolviendo el lecho.
"""

import os
import sys
import subprocess

# Auto-detección y re-ejecución con el entorno virtual .venv
root_dir = os.path.dirname(os.path.abspath(__file__))
venv_python = os.path.join(root_dir, ".venv", "Scripts", "python.exe")

if os.path.exists(venv_python) and os.path.normcase(sys.executable) != os.path.normcase(venv_python):
    result = subprocess.run([venv_python] + sys.argv)
    sys.exit(result.returncode)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

R_UNIVERSAL = 8.31446  # J/(mol·K)

def menu_configuracion():
    print("==================================================================")
    print("   CONFIGURACIÓN INTERACTIVA: HORNO + REACTOR DE ACERO + LECHO    ")
    print("==================================================================")
    print("Presiona ENTER para usar el valor sugerido entre corchetes [ ].\n")

    def pedir_valor(mensaje, default, tipo=float):
        resp = input(f"{mensaje} [{default}]: ").strip()
        if not resp:
            return default
        try:
            return tipo(resp)
        except ValueError:
            print(f"  -> Usando valor por defecto: {default}")
            return default

    ramp_rate = pedir_valor("1. Velocidad de rampa de calentamiento (°C/min)", 10.0)
    T_final   = pedir_valor("2. Temperatura máxima de pirólisis (°C)", 500.0)
    soak_time = pedir_valor("3. Tiempo de sostenimiento isotérmico / Soak (min)", 30.0)
    tau_stone = pedir_valor("4. Inercia de la piedra refractaria tau (s)", 200.0)
    R_bed_mm  = pedir_valor("5. Radio interior de la canastilla de café (mm)", 25.0)
    espesor_acero_mm = pedir_valor("6. Espesor del tubo de acero inoxidable (mm)", 2.5)
    L_bed_cm  = pedir_valor("7. Altura de la canastilla de café (cm)", 20.0)
    humedad   = pedir_valor("8. Humedad del café acondicionado (% peso)", 8.0)

    return {
        'ramp_rate': ramp_rate,
        'T_final': T_final,
        'soak_time': soak_time,
        'tau_stone': tau_stone,
        'R_int': R_bed_mm / 1000.0,
        'e_inox': espesor_acero_mm / 1000.0,
        'L': L_bed_cm / 100.0,
        'w_H2O': humedad / 100.0
    }

def ejecutar_simulacion(cfg):
    R_int = cfg['R_int']
    e_inox = cfg['e_inox']
    R_ext = R_int + e_inox
    L = cfg['L']

    # Discretización radial multicapa:
    # 0 a R_int: Lecho de café (Nr_bed nodos)
    # R_int a R_ext: Pared de acero inox (Nr_inox nodos)
    Nr_bed = 14
    Nr_inox = 4
    Nr_total = Nr_bed + Nr_inox
    Nz = 20

    r_bed = np.linspace(0, R_int, Nr_bed)
    r_inox = np.linspace(R_int + (e_inox / Nr_inox), R_ext, Nr_inox)
    r_grid = np.concatenate([r_bed, r_inox])
    z_grid = np.linspace(0, L, Nz)
    dr_bed = r_bed[1] - r_bed[0]
    dr_inox = r_inox[1] - r_inox[0] if Nr_inox > 1 else e_inox
    dz = L / (Nz - 1)

    # Volúmenes diferenciales del café
    dV_bed = np.zeros((Nr_bed, Nz))
    for i in range(Nr_bed):
        r_i = r_bed[i]
        if i == 0:
            area = np.pi * (dr_bed / 2.0)**2
        elif i == Nr_bed - 1:
            area = np.pi * (R_int**2 - (R_int - dr_bed / 2.0)**2)
        else:
            area = 2.0 * np.pi * r_i * dr_bed
        dV_bed[i, :] = area * dz

    # Propiedades del café y del acero
    rho_bulk = 380.0
    w_H2O = cfg['w_H2O']
    w_ASH = 0.0231 * (1.0 - w_H2O)
    rho_water_0 = rho_bulk * w_H2O
    rho_char_0  = rho_bulk * w_ASH
    rho_bio_0   = rho_bulk - rho_water_0 - rho_char_0

    rho_inox = 7900.0
    cp_inox  = 500.0
    k_inox   = 16.3

    T_ini_K = 25.0 + 273.15
    T_fin_K = cfg['T_final'] + 273.15
    beta_K_s = cfg['ramp_rate'] / 60.0
    t_ramp_end = (T_fin_K - T_ini_K) / beta_K_s
    t_final = t_ramp_end + (cfg['soak_time'] * 60.0)

    # Cinética de pirólisis calibrada
    A_char = 4.2e5;        Ea_char = 88.0e3;     dH_char = 120.0e3
    A_vol  = 2.8e10;       Ea_vol  = 155.0e3;    dH_vol  = -210.0e3
    A_gnc  = 8.0e6;        Ea_gnc  = 105.0e3;    dH_gnc  = -80.0e3
    A_evap = 5.0e5;        Ea_evap = 41.8e3;     dH_evap = -2.26e6

    factor_sec_char = np.clip(0.18 / np.sqrt(max(cfg['ramp_rate'] / 10.0, 0.2)), 0.04, 0.25)

    # Matrices de estado
    T       = np.full((Nr_total, Nz), T_ini_K)
    biomass = np.full((Nr_bed, Nz), rho_bio_0)
    char    = np.full((Nr_bed, Nz), rho_char_0)
    water   = np.full((Nr_bed, Nz), rho_water_0)

    M_init_total = np.sum((biomass + char + water) * dV_bed)

    time_hist = []
    T_center_hist = []
    T_mid_hist = []
    T_bed_surf_hist = []
    T_inox_ext_hist = []
    T_mean_bed_hist = []
    T_furnace_hist = []
    T_set_hist = []

    pct_biomass_hist = []
    pct_char_hist = []
    pct_volatiles_hist = []
    pct_gnc_hist = []
    pct_water_hist = []
    residual_mass_hist = []

    snapshots_T = []
    snapshot_times = []
    save_interval = max(5.0, t_final / 120.0)
    last_save = -save_interval

    cum_vol = 0.0
    cum_gnc = 0.0
    cum_vap = 0.0

    t = 0.0
    T_furnace = T_ini_K
    tau_stone = cfg['tau_stone']

    print(f"\nIniciando simulación con reactor de acero inox (R_ext = {R_ext*1000:.1f} mm, e = {e_inox*1000:.1f} mm)...")
    while t <= t_final:
        # Horno con inercia de piedra refractaria
        T_set = T_ini_K + beta_K_s * t if t <= t_ramp_end else T_fin_K
        T_furnace += ((T_set - T_furnace) / max(tau_stone, 10.0)) * 0.4

        alpha_bed = 0.12 / (rho_bulk * 1500.0)
        cfl_dt = 0.35 * (dr_bed**2) / alpha_bed
        dt = min(0.3, cfl_dt, t_final - t if t_final > t else 0.3)
        if dt <= 0:
            break

        # Propiedades del lecho
        conv = np.clip(1.0 - (biomass + water) / (rho_bio_0 + rho_water_0), 0.0, 1.0)
        k_bed = (1.0 - conv) * (0.12 + 1.5e-4 * (T[:Nr_bed, :] - 298.15)) + conv * (0.08 + 2.0e-4 * (T[:Nr_bed, :] - 298.15))
        cp_bed = (1.0 - conv) * (1250.0 + 1.8 * (T[:Nr_bed, :] - 273.15)) + conv * (720.0 + 2.2 * (T[:Nr_bed, :] - 273.15))
        rho_bed = np.maximum(biomass + char + water, 50.0)

        # Cinética de pirólisis en el lecho
        T_bed_clip = np.maximum(T[:Nr_bed, :], 273.15)
        RT_bed = R_UNIVERSAL * T_bed_clip
        k_c = A_char * np.exp(-Ea_char / RT_bed)
        k_v = A_vol  * np.exp(-Ea_vol  / RT_bed)
        k_g = A_gnc  * np.exp(-Ea_gnc  / RT_bed)
        k_w = np.where(T_bed_clip > 350.0, A_evap * np.exp(-Ea_evap / RT_bed), 0.0)

        r_c_prim = k_c * np.maximum(biomass, 0.0)
        r_v_prim = k_v * np.maximum(biomass, 0.0)
        r_g_prim = k_g * np.maximum(biomass, 0.0)
        r_w      = k_w * np.maximum(water, 0.0)

        r_c_sec = factor_sec_char * r_v_prim
        r_v_net = r_v_prim - r_c_sec

        d_bio  = -(r_c_prim + r_v_prim + r_g_prim)
        d_char = r_c_prim + r_c_sec
        d_wat  = -r_w
        q_rxn_bed = r_w * dH_evap + r_c_prim * dH_char + r_v_prim * dH_vol + r_g_prim * dH_gnc

        # Conducción en todo el dominio (Lecho + Acero)
        dT_dt = np.zeros_like(T)
        dz2 = dz**2

        # A. Centro del lecho r = 0 (L'Hôpital)
        term_rad_center = 4.0 * k_bed[0, :] * (T[1, :] - T[0, :]) / (dr_bed**2)
        d2T_dz2_center = np.zeros(Nz)
        d2T_dz2_center[1:-1] = (T[0, 2:] - 2*T[0, 1:-1] + T[0, :-2]) / dz2
        dT_dt[0, :] = (term_rad_center + k_bed[0, :] * d2T_dz2_center + q_rxn_bed[0, :]) / (rho_bed[0, :] * cp_bed[0, :])

        # B. Nodos internos del lecho (1 <= i <= Nr_bed - 2)
        for i in range(1, Nr_bed - 1):
            d2T_dr2 = (T[i+1, :] - 2*T[i, :] + T[i-1, :]) / (dr_bed**2)
            dT_dr   = (T[i+1, :] - T[i-1, :]) / (2.0 * dr_bed)
            d2T_dz2 = np.zeros(Nz)
            d2T_dz2[1:-1] = (T[i, 2:] - 2*T[i, 1:-1] + T[i, :-2]) / dz2
            term_rad = k_bed[i, :] * (d2T_dr2 + (1.0 / r_bed[i]) * dT_dr)
            dT_dt[i, :] = (term_rad + k_bed[i, :] * d2T_dz2 + q_rxn_bed[i, :]) / (rho_bed[i, :] * cp_bed[i, :])

        # C. Interfaz Lecho - Acero Inoxidable (i = Nr_bed - 1)
        idx_int = Nr_bed - 1
        q_transf_int = (T[idx_int + 1, :] - T[idx_int, :]) / ((dr_bed / (2 * k_bed[-1, :])) + (dr_inox / (2 * k_inox)) + 0.003)
        d2T_dr2_int = 2.0 * (T[idx_int - 1, :] - T[idx_int, :] + (dr_bed * q_transf_int / k_bed[-1, :])) / (dr_bed**2)
        term_rad_int = k_bed[-1, :] * (d2T_dr2_int + (1.0 / R_int) * (q_transf_int / k_bed[-1, :]))
        d2T_dz2_int = np.zeros(Nz)
        d2T_dz2_int[1:-1] = (T[idx_int, 2:] - 2*T[idx_int, 1:-1] + T[idx_int, :-2]) / dz2
        dT_dt[idx_int, :] = (term_rad_int + k_bed[-1, :] * d2T_dz2_int + q_rxn_bed[-1, :]) / (rho_bed[-1, :] * cp_bed[-1, :])

        # D. Nodos dentro del espesor del Acero Inoxidable
        for i in range(Nr_bed, Nr_total - 1):
            r_curr = r_grid[i]
            d2T_dr2 = (T[i+1, :] - 2*T[i, :] + T[i-1, :]) / (dr_inox**2)
            dT_dr   = (T[i+1, :] - T[i-1, :]) / (2.0 * dr_inox)
            d2T_dz2 = np.zeros(Nz)
            d2T_dz2[1:-1] = (T[i, 2:] - 2*T[i, 1:-1] + T[i, :-2]) / dz2
            term_rad = k_inox * (d2T_dr2 + (1.0 / r_curr) * dT_dr)
            dT_dt[i, :] = (term_rad + k_inox * d2T_dz2) / (rho_inox * cp_inox)

        # E. Cara externa del reactor de acero r = R_ext (Robin con horno)
        h_comb = 90.0
        q_horno = h_comb * (T_furnace - T[-1, :])
        d2T_dr2_ext = 2.0 * (T[-2, :] - T[-1, :] + (dr_inox * q_horno / k_inox)) / (dr_inox**2)
        term_rad_ext = k_inox * (d2T_dr2_ext + (1.0 / R_ext) * (q_horno / k_inox))
        d2T_dz2_ext = np.zeros(Nz)
        d2T_dz2_ext[1:-1] = (T[-1, 2:] - 2*T[-1, 1:-1] + T[-1, :-2]) / dz2
        dT_dt[-1, :] = (term_rad_ext + k_inox * d2T_dz2_ext) / (rho_inox * cp_inox)

        # Integración
        T += dT_dt * dt
        biomass = np.maximum(biomass + d_bio * dt, 0.0)
        char    = np.maximum(char + d_char * dt, 0.0)
        water   = np.maximum(water + d_wat * dt, 0.0)

        cum_vol += np.sum(r_v_net * dV_bed) * dt
        cum_gnc += np.sum(r_g_prim * dV_bed) * dt
        cum_vap += np.sum(r_w * dV_bed) * dt

        if (t - last_save) >= save_interval or t >= t_final:
            last_save = t
            time_hist.append(t / 60.0)
            mid_z = Nz // 2
            T_center_hist.append(T[0, mid_z] - 273.15)
            T_mid_hist.append(T[Nr_bed // 2, mid_z] - 273.15)
            T_bed_surf_hist.append(T[Nr_bed - 1, mid_z] - 273.15)
            T_inox_ext_hist.append(T[-1, mid_z] - 273.15)
            T_mean_bed_hist.append(np.mean(T[:Nr_bed, :]) - 273.15)
            T_furnace_hist.append(T_furnace - 273.15)
            T_set_hist.append(T_set - 273.15)

            M_bio_curr  = np.sum(biomass * dV_bed)
            M_char_curr = np.sum(char * dV_bed)
            pct_biomass_hist.append((M_bio_curr / M_init_total) * 100.0)
            pct_char_hist.append((M_char_curr / M_init_total) * 100.0)
            pct_volatiles_hist.append((cum_vol / M_init_total) * 100.0)
            pct_gnc_hist.append((cum_gnc / M_init_total) * 100.0)
            pct_water_hist.append((cum_vap / M_init_total) * 100.0)
            residual_mass_hist.append(((M_bio_curr + M_char_curr + np.sum(water * dV_bed)) / M_init_total) * 100.0)

            snapshots_T.append(np.copy(T) - 273.15)
            snapshot_times.append(t / 60.0)

        t += dt

    dt_arr = np.diff(time_hist)
    dW_arr = np.diff(residual_mass_hist)
    dtg = np.zeros_like(time_hist)
    dtg[1:] = -dW_arr / np.maximum(dt_arr, 1e-4)

    return {
        'time': np.array(time_hist),
        'T_center': np.array(T_center_hist),
        'T_mid': np.array(T_mid_hist),
        'T_bed_surf': np.array(T_bed_surf_hist),
        'T_inox_ext': np.array(T_inox_ext_hist),
        'T_mean': np.array(T_mean_bed_hist),
        'T_furnace': np.array(T_furnace_hist),
        'T_set': np.array(T_set_hist),
        'pct_biomass': np.array(pct_biomass_hist),
        'pct_char': np.array(pct_char_hist),
        'pct_vol': np.array(pct_volatiles_hist),
        'pct_gnc': np.array(pct_gnc_hist),
        'pct_water': np.array(pct_water_hist),
        'residual_mass': np.array(residual_mass_hist),
        'dtg': dtg,
        'snapshots_T': snapshots_T,
        'snapshot_times': snapshot_times,
        'r_grid_mm': r_grid * 1000.0,
        'z_grid_cm': z_grid * 100.0,
        'R_int_mm': R_int * 1000.0,
        'R_ext_mm': R_ext * 1000.0,
        'Nr_bed': Nr_bed,
        'cfg': cfg
    }

def graficar_resultados(res):
    t = res['time']
    T_mean = res['T_mean']

    # =========================================================================
    # VENTANA 1: PANEL CIENTÍFICO DE 6 FIGURAS
    # =========================================================================
    fig1 = plt.figure(figsize=(17, 10.5))
    fig1.canvas.manager.set_window_title("Análisis Completo: Horno + Reactor de Acero + Lecho")

    # 1. Temperaturas vs Tiempo
    ax1 = plt.subplot(2, 3, 1)
    ax1.plot(t, res['T_set'], 'k--', label="Setpoint (PID)")
    ax1.plot(t, res['T_furnace'], 'm-', label="Horno (Piedra)")
    ax1.plot(t, res['T_inox_ext'], color='darkorange', linewidth=2.5, label=f"Exterior Tubo Inox (r={res['R_ext_mm']:.1f}mm)")
    ax1.plot(t, res['T_bed_surf'], 'r-', label=f"Superficie Café (r={res['R_int_mm']:.1f}mm)")
    ax1.plot(t, res['T_mid'], 'g-.', label="Café radio medio")
    ax1.plot(t, res['T_center'], 'b-', label="Centro Café (r=0, L'Hôpital)")
    ax1.set_xlabel("Tiempo (min)")
    ax1.set_ylabel("Temperatura (°C)")
    ax1.set_title("Evolución Térmica Multicapa", fontweight='bold')
    ax1.legend(loc='lower right', fontsize=8)
    ax1.grid(True, linestyle=':', alpha=0.6)

    # 2. Especies vs Tiempo
    ax2 = plt.subplot(2, 3, 2)
    ax2.plot(t, res['pct_biomass'], 'brown', linewidth=2, label="Biomasa virgen")
    ax2.plot(t, res['pct_char'], 'black', linewidth=2, label="Biochar")
    ax2.plot(t, res['pct_vol'], 'orange', linewidth=2, label="Biovolátiles / Bio-oil")
    ax2.plot(t, res['pct_gnc'], 'blue', linewidth=1.8, label="Gases No Condensables (GNC)")
    ax2.plot(t, res['pct_water'], 'teal', linestyle='--', label="Vapor H2O")
    ax2.set_xlabel("Tiempo (min)")
    ax2.set_ylabel("Fracción en masa (% wt)")
    ax2.set_title("Evolución de Especies vs Tiempo", fontweight='bold')
    ax2.legend(loc='center left', fontsize=8)
    ax2.grid(True, linestyle=':', alpha=0.6)

    # 3. Especies vs Temperatura Media
    ax3 = plt.subplot(2, 3, 3)
    ax3.plot(T_mean, res['pct_biomass'], 'brown', linewidth=2, label="Biomasa virgen")
    ax3.plot(T_mean, res['pct_char'], 'black', linewidth=2, label="Biochar")
    ax3.plot(T_mean, res['pct_vol'], 'orange', linewidth=2, label="Biovolátiles / Bio-oil")
    ax3.plot(T_mean, res['pct_gnc'], 'blue', linewidth=1.8, label="Gases No Condensables (GNC)")
    ax3.set_xlabel("Temperatura media del lecho (°C)")
    ax3.set_ylabel("Fracción en masa (% wt)")
    ax3.set_title("Especies vs Temperatura del Lecho", fontweight='bold')
    ax3.legend(loc='center left', fontsize=8)
    ax3.grid(True, linestyle=':', alpha=0.6)

    # 4. TGA y DTG
    ax4 = plt.subplot(2, 3, 4)
    ax4_dtg = ax4.twinx()
    l1 = ax4.plot(T_mean, res['residual_mass'], 'k-', linewidth=2, label="TGA: Masa residual (% wt)")
    l2 = ax4_dtg.plot(T_mean, res['dtg'], 'r--', linewidth=1.5, label="DTG: -dW/dt (%/min)")
    ax4.set_xlabel("Temperatura media (°C)")
    ax4.set_ylabel("Masa residual (% peso)", color='black')
    ax4_dtg.set_ylabel("DTG: Pérdida (%/min)", color='red')
    ax4.set_title("Termogravimetría (TGA y DTG)", fontweight='bold')
    ax4.set_ylim([0, 105])
    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax4.legend(lines, labels, loc='center right', fontsize=8)
    ax4.grid(True, linestyle=':', alpha=0.6)

    # 5. Gradientes Térmicos a través del Acero y el Café
    ax5 = plt.subplot(2, 3, 5)
    delta_T_acero = res['T_inox_ext'] - res['T_bed_surf']
    delta_T_cafe  = res['T_bed_surf'] - res['T_center']
    ax5.plot(t, delta_T_acero, color='darkorange', linewidth=2, label="ΔT a través del Tubo de Acero")
    ax5.plot(t, delta_T_cafe, color='purple', linewidth=2, label="ΔT a través del Café (Pared-Centro)")
    ax5.set_xlabel("Tiempo (min)")
    ax5.set_ylabel("Diferencia de Temperatura (°C)")
    ax5.set_title("Gradientes Térmicos en las Interfaces", fontweight='bold')
    ax5.legend(loc='upper right', fontsize=8)
    ax5.grid(True, linestyle=':', alpha=0.6)

    # 6. Diagrama de barras de rendimientos finales
    ax6 = plt.subplot(2, 3, 6)
    categorias = ['Biochar\n(Sólido)', 'Biovolátiles\n(Bio-oil)', 'GNC\n(Gases)', 'Vapor\n(Agua)']
    valores = [res['pct_char'][-1], res['pct_vol'][-1], res['pct_gnc'][-1], res['pct_water'][-1]]
    colores = ['#3E2723', '#BF360C', '#1A237E', '#00695C']
    barras = ax6.bar(categorias, valores, color=colores, edgecolor='black', width=0.55)
    ax6.set_ylabel("Rendimiento final (% wt)")
    ax6.set_title(f"Rendimientos Finales (Rampa {res['cfg']['ramp_rate']} °C/min)", fontweight='bold')
    ax6.set_ylim([0, max(valores) * 1.25])
    for b in barras:
        h = b.get_height()
        ax6.annotate(f"{h:.1f}%", xy=(b.get_x() + b.get_width()/2, h), xytext=(0, 3),
                     textcoords="offset points", ha='center', va='bottom', fontweight='bold')
    ax6.grid(axis='y', linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig("analisis_completo_pirolisis.png", dpi=300)
    plt.show(block=False)

    # =========================================================================
    # VENTANA 2: DOS MAPAS TÉRMICOS 2D INTERACTIVOS CON SLIDER
    # =========================================================================
    fig_map = plt.figure(figsize=(15, 8))
    fig_map.canvas.manager.set_window_title("Mapas Térmicos 2D Interactivos (Tubo de Acero + Lecho)")
    plt.subplots_adjust(bottom=0.18, wspace=0.30)

    # Malla para el mapa longitudinal
    R_mesh, Z_mesh = np.meshgrid(res['r_grid_mm'], res['z_grid_cm'], indexing='ij')
    vmin = min([np.min(snap) for snap in res['snapshots_T']])
    vmax = max([np.max(snap) for snap in res['snapshots_T']])

    # Malla polar para el corte circular transversal (plano x, y)
    theta = np.linspace(0, 2*np.pi, 60)
    R_polar_mesh, Theta_mesh = np.meshgrid(res['r_grid_mm'], theta, indexing='ij')
    X_polar = R_polar_mesh * np.cos(Theta_mesh)
    Y_polar = R_polar_mesh * np.sin(Theta_mesh)

    # Subplot 1: Mapa Longitudinal (r vs z)
    ax_map1 = fig_map.add_subplot(1, 2, 1)
    # Subplot 2: Corte Transversal Circular (x vs y)
    ax_map2 = fig_map.add_subplot(1, 2, 2)

    mid_z = len(res['z_grid_cm']) // 2

    def dibujar_mapas(idx):
        ax_map1.cla()
        ax_map2.cla()
        snap = res['snapshots_T'][idx]

        # 1. Mapa Longitudinal
        cs1 = ax_map1.contourf(R_mesh, Z_mesh, snap, levels=np.linspace(vmin, vmax, 40), cmap='inferno')
        ax_map1.axvline(x=res['R_int_mm'], color='cyan', linestyle='--', linewidth=2, label="Límite: Acero / Café")
        ax_map1.set_xlabel("Radio r (mm)", fontweight='bold')
        ax_map1.set_ylabel("Altura axial z (cm)", fontweight='bold')
        ax_map1.set_title(f"Vista Longitudinal 2D (r vs z)\n[0 a {res['R_int_mm']:.0f} mm: Café | {res['R_int_mm']:.0f} a {res['R_ext_mm']:.1f} mm: Acero Inox]", fontsize=10, fontweight='bold')
        ax_map1.legend(loc='lower left', fontsize=8)

        # 2. Corte Transversal Circular
        T_radial_profile = snap[:, mid_z]
        T_polar_matrix = np.tile(T_radial_profile[:, np.newaxis], (1, len(theta)))
        cs2 = ax_map2.contourf(X_polar, Y_polar, T_polar_matrix, levels=np.linspace(vmin, vmax, 40), cmap='inferno')
        # Círculo de la interfaz
        circ_int = plt.Circle((0, 0), res['R_int_mm'], color='cyan', fill=False, linestyle='--', linewidth=2, label="Canastilla interior")
        circ_ext = plt.Circle((0, 0), res['R_ext_mm'], color='white', fill=False, linestyle='-', linewidth=1.5, label="Exterior tubo inox")
        ax_map2.add_patch(circ_int)
        ax_map2.add_patch(circ_ext)
        ax_map2.set_aspect('equal', 'box')
        ax_map2.set_xlabel("Coordenada X (mm)", fontweight='bold')
        ax_map2.set_ylabel("Coordenada Y (mm)", fontweight='bold')
        ax_map2.set_title(f"Corte Circular Transversal a media altura (z = {res['z_grid_cm'][mid_z]:.1f} cm)\n[Anillo Exterior = Reactor de Acero Inox]", fontsize=10, fontweight='bold')
        ax_map2.legend(loc='lower right', fontsize=8)

        return cs1

    cs1 = dibujar_mapas(0)

    # Barra de color común
    cbar_ax = fig_map.add_axes([0.25, 0.92, 0.50, 0.025])
    cbar = fig_map.colorbar(cs1, cax=cbar_ax, orientation='horizontal')
    cbar.set_label("Temperatura (°C)", fontweight='bold')

    # Slider interactivo
    ax_slider = fig_map.add_axes([0.25, 0.05, 0.50, 0.03], facecolor='lightgoldenrodyellow')
    slider_tiempo = Slider(ax_slider, 'Tiempo (min)', 0, len(res['snapshot_times']) - 1, valinit=0, valstep=1)

    def actualizar(val):
        idx = int(slider_tiempo.val)
        dibujar_mapas(idx)
        slider_tiempo.valtext.set_text(f"{res['snapshot_times'][idx]:.1f} min")
        fig_map.canvas.draw_idle()

    slider_tiempo.on_changed(actualizar)
    plt.suptitle(f"Simulación Térmica Multicapa | Tiempo: {res['snapshot_times'][0]:.1f} min", y=0.99, fontsize=12, fontweight='bold')
    plt.show()

if __name__ == "__main__":
    cfg = menu_configuracion()
    resultados = ejecutar_simulacion(cfg)
    print("\n------------------------------------------------------------------")
    print(f"RENDIMIENTOS FINALES ESTIMADOS (Rampa = {cfg['ramp_rate']} °C/min):")
    print(f"  • Biochar (sólido carbonoso):              {resultados['pct_char'][-1]:5.1f} %")
    print(f"  • Biovolátiles / Bio-oil (condensables):   {resultados['pct_vol'][-1]:5.1f} %")
    print(f"  • Gases No Condensables (GNC):             {resultados['pct_gnc'][-1]:5.1f} %")
    print(f"  • Vapor de agua (humedad):                 {resultados['pct_water'][-1]:5.1f} %")
    total_balance = resultados['pct_char'][-1] + resultados['pct_vol'][-1] + resultados['pct_gnc'][-1] + resultados['pct_water'][-1]
    print(f"  TOTAL BALANCE:                             {total_balance:5.1f} %")
    print("------------------------------------------------------------------")
    graficar_resultados(resultados)