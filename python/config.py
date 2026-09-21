"""
config.py
Configuración de parámetros geométricos, termofísicos, cinéticos y operativos
para el reactor vertical de pirólisis de residuos de café (SCG).
Basado en datos experimentales de laboratorio (GIEM - UdeA / EIA) y espectroscopía FTIR.
"""

import numpy as np

# ==============================================================================
# 1. PARÁMETROS GEOMÉTRICOS DEL REACTOR Y DE LA CANASTILLA
# ==============================================================================
class GeometryConfig:
    def __init__(
        self,
        R_inner=0.025,       # Radio interno del lecho de biomasa (m) [50 mm diámetro]
        R_outer_tube=0.027,  # Radio externo del tubo de acero inox (m) [espesor pared = 2 mm]
        L_bed=0.20,          # Longitud/altura de la canastilla de café (m) [20 cm]
        Nr=15,               # Número de nodos en dirección radial (r)
        Nz=20                # Número de nodos en dirección axial (z)
    ):
        self.R = R_inner
        self.R_tube = R_outer_tube
        self.L = L_bed
        self.Nr = Nr
        self.Nz = Nz
        self.dr = self.R / (self.Nr - 1)
        self.dz = self.L / (self.Nz - 1)
        self.r_grid = np.linspace(0.0, self.R, self.Nr)
        self.z_grid = np.linspace(0.0, self.L, self.Nz)


# ==============================================================================
# 2. CARACTERIZACIÓN DEL RESIDUO DE CAFÉ (GIEM + FTIR)
# ==============================================================================
class CoffeeBiomassConfig:
    def __init__(
        self,
        # Densidad aparente empacada en la canastilla
        rho_bulk=380.0,      # kg/m^3 (lecho poroso empacado de café seco/molido)
        # Humedad del café al ingresar al reactor (base húmeda, % wt)
        # El café fresco tiene 61.2% (reporte GIEM), pero antes de pirolizar se seca
        # a estufa o al sol típicamente a un 5% - 10%.
        moisture_fraction=0.08,
        # Composición lignocelulósica en base seca (% wt), coherente con FTIR y GIEM:
        # FTIR revela picos intensos alifáticos de aceites de café a 2924 y 2853 cm^-1,
        # carbonilo éster a 1743 cm^-1, polisacáridos (1030-1160 cm^-1) y aromáticos de lignina (1500-1650 cm^-1).
        frac_cellulose=0.12,      # Celulosa (~12%)
        frac_hemicellulose=0.35,  # Hemicelulosa (galactomananos) (~35%)
        frac_lignin=0.28,         # Lignina total (~28%)
        frac_lipids=0.15,         # Aceites/lípidos de café (~15%)
        frac_ash=0.0231           # Cenizas base seca (2.31%, Reporte GIEM)
    ):
        self.rho_bulk = rho_bulk
        self.w_H2O = moisture_fraction
        self.w_CELL = frac_cellulose * (1.0 - self.w_H2O)
        self.w_HCE = frac_hemicellulose * (1.0 - self.w_H2O)
        self.w_LIG = frac_lignin * (1.0 - self.w_H2O)
        self.w_LIP = frac_lipids * (1.0 - self.w_H2O)
        self.w_ASH = frac_ash * (1.0 - self.w_H2O)
        
        # Porosidad inicial del lecho de partículas
        self.eps_solid_0 = 0.45
        self.eps_gas_0 = 1.0 - self.eps_solid_0
        
        # Densidades reales intrínsecas (kg/m^3)
        self.rho_cell = 1500.0
        self.rho_hce = 1400.0
        self.rho_lig = 1350.0
        self.rho_char = 1300.0
        self.rho_water = 1000.0
        self.rho_gas = 0.8
        self.rho_tar = 1.1


# ==============================================================================
# 3. PROPIEDADES TERMOFÍSICAS Y TRANSPORTE
# ==============================================================================
class ThermalPropertiesConfig:
    def __init__(self):
        # Conductividad térmica efectiva de la biomasa porosa (W/(m·K))
        # k(T) = k0 + a*T
        self.k_biomass_0 = 0.12
        self.k_biomass_slope = 1.5e-4
        self.k_char_0 = 0.08
        self.k_char_slope = 2.0e-4
        self.k_gas = 0.035
        self.k_steel = 16.3   # Tubo de acero inoxidable AISI 304 (W/(m·K))
        
        # Capacidad calorífica específica Cp (J/(kg·K))
        # Funciones dependientes de T
        self.cp_water_liquid = 4184.0
        self.cp_steam = 2050.0
        self.cp_steel = 500.0
        
        # Coeficiente combinado de convección y radiación horno-tubo (W/(m^2·K))
        self.h_furnace_tube = 85.0
        
        # Permeabilidad y viscosidad de gases pirolíticos (flujo de Darcy)
        self.permeability_bed = 5e-11  # m^2 (lecho granular poroso)
        self.visc_gas = 3.2e-5         # Pa·s


# ==============================================================================
# 4. HORNO ELÉCTRICO Y CONTROLADOR PID (MaxThermo MC-5438)
# ==============================================================================
class FurnaceConfig:
    def __init__(
        self,
        power_kW=2.5,             # Potencia nominal de resistencias calefactoras (kW)
        stone_thermal_inertia=240.0, # Constante de tiempo tau (s) de la piedra refractaria
        h_loss_ambient=10.0,      # Pérdidas al ambiente exterior (W/(m^2·K))
        T_ambient=298.15          # Temperatura ambiente (25 °C = 298.15 K)
    ):
        self.power_max = power_kW * 1000.0
        self.tau_stone = stone_thermal_inertia
        self.h_loss = h_loss_ambient
        self.T_amb = T_ambient


# ==============================================================================
# 5. PROGRAMACIÓN DE RAMPAS TÉRMICAS
# ==============================================================================
class HeatingRampConfig:
    """
    Define el programa de temperatura del horno:
    - Rampas de subida en °C/min (ej. 5 °C/min, 10 °C/min, 20 °C/min)
    - Mesetas isotérmicas (soak time) a temperatura final de pirólisis (ej. 450 - 600 °C)
    - Enfriamiento opcional
    """
    def __init__(
        self,
        T_initial_C=25.0,
        ramp_rate_C_per_min=10.0,
        T_final_C=500.0,
        soak_time_min=45.0
    ):
        self.T_init_K = T_initial_C + 273.15
        self.T_final_K = T_final_C + 273.15
        self.ramp_rate_K_s = ramp_rate_C_per_min / 60.0 # K/s
        
        # Tiempo requerido para alcanzar T_final
        self.t_ramp_end = (self.T_final_K - self.T_init_K) / self.ramp_rate_K_s
        self.t_soak_end = self.t_ramp_end + (soak_time_min * 60.0)
        self.total_time = self.t_soak_end

    def get_setpoint_T(self, t):
        """Devuelve la temperatura objetivo (K) del controlador en el instante t (segundos)"""
        if t <= 0:
            return self.T_init_K
        elif t <= self.t_ramp_end:
            return self.T_init_K + self.ramp_rate_K_s * t
        else:
            return self.T_final_K
