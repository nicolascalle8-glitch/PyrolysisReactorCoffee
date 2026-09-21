"""
kinetics.py
Mecanismo cinético de pirólisis multicomponente para residuos de café (SCG).
Modela la degradación térmica secuencial y paralela de:
1. Humedad libre y ligada (secado / evaporación)
2. Hemicelulosa (200 - 320 °C)
3. Celulosa (300 - 400 °C)
4. Lignina (250 - 650 °C)
5. Aceites/lípidos de café (350 - 480 °C, identificados por FTIR en 2924 y 2853 cm^-1)
6. Reacciones secundarias de craqueo de alquitranes (tar cracking -> gas + biochar secundario)
"""

import numpy as np

R_UNIVERSAL = 8.31446  # J/(mol·K)

class PyrolysisKinetics:
    def __init__(self):
        # ----------------------------------------------------------------------
        # Factores pre-exponenciales A (s^-1) y Energías de activación Ea (J/mol)
        # Adaptados de modelos cinéticos de biomasa (Ranzi, Di Blasi, Neves) y SCG
        # ----------------------------------------------------------------------
        # Evaporación de agua
        self.A_evap = 5.0e5
        self.Ea_evap = 41.8e3
        self.dH_evap = -2.26e6  # J/kg (endotérmico)
        
        # Hemicelulosa (HCE)
        self.A_hce = 1.2e10
        self.Ea_hce = 125.0e3
        self.dH_hce = -65.0e3   # J/kg (ligeramente endotérmico)
        # Rendimientos estequiométricos primarios HCE
        self.y_char_hce = 0.28
        self.y_tar_hce  = 0.38
        self.y_gas_hce  = 0.34
        
        # Celulosa (CELL)
        self.A_cell = 2.5e13
        self.Ea_cell = 185.0e3
        self.dH_cell = -150.0e3  # J/kg (endotérmico)
        # Rendimientos primarios CELL
        self.y_char_cell = 0.12
        self.y_tar_cell  = 0.65
        self.y_gas_cell  = 0.23
        
        # Lignina (LIG)
        self.A_lig = 8.0e7
        self.Ea_lig = 98.0e3
        self.dH_lig = 120.0e3   # J/kg (ligeramente exotérmico debido a policondensación)
        # Rendimientos primarios LIG
        self.y_char_lig = 0.52
        self.y_tar_lig  = 0.28
        self.y_gas_lig  = 0.20
        
        # Lípidos / Aceites de café (LIP)
        self.A_lip = 2.0e9
        self.Ea_lip = 120.0e3
        self.dH_lip = -180.0e3
        # Rendimientos primarios LIP (alta volatilización a bio-oil)
        self.y_char_lip = 0.05
        self.y_tar_lip  = 0.82
        self.y_gas_lip  = 0.13
        
        # Craqueo secundario de volátiles/alquitranes (Tar cracking)
        # Tar -> Gas (reacción homogénea)
        self.A_crk_gas = 4.2e6
        self.Ea_crk_gas = 108.0e3
        self.dH_crk_gas = -42.0e3
        
        # Tar -> Char secundario + Gas (reacción heterogénea catalizada por char)
        self.A_crk_char = 1.0e5
        self.Ea_crk_char = 95.0e3
        self.dH_crk_char = 60.0e3

    def compute_rates(self, T, rho_w, rho_hce, rho_cell, rho_lig, rho_lip, rho_tar, rho_char):
        """
        Calcula las velocidades de reacción de cada componente y la tasa neta
        de generación de calor por pirólisis (q_dot_rxn).
        
        T: Temperatura local (K)
        rho_*: Concentraciones másicas locales en kg/m^3 de lecho
        """
        T_clipped = np.maximum(T, 250.0) # Prevenir overflow numérico a bajas T
        RT = R_UNIVERSAL * T_clipped
        
        # Constantes de Arrhenius
        k_evap = np.where(T_clipped > 350.0, self.A_evap * np.exp(-self.Ea_evap / RT), 0.0)
        k_hce  = self.A_hce  * np.exp(-self.Ea_hce  / RT)
        k_cell = self.A_cell * np.exp(-self.Ea_cell / RT)
        k_lig  = self.A_lig  * np.exp(-self.Ea_lig  / RT)
        k_lip  = self.A_lip  * np.exp(-self.Ea_lip  / RT)
        
        k_crk_g = self.A_crk_gas  * np.exp(-self.Ea_crk_gas  / RT)
        k_crk_c = self.A_crk_char * np.exp(-self.Ea_crk_char / RT)
        
        # Tasas de consumo de sólidos y líquidos (kg/(m^3·s))
        r_evap = k_evap * np.maximum(rho_w, 0.0)
        r_hce  = k_hce  * np.maximum(rho_hce, 0.0)
        r_cell = k_cell * np.maximum(rho_cell, 0.0)
        r_lig  = k_lig  * np.maximum(rho_lig, 0.0)
        r_lip  = k_lip  * np.maximum(rho_lip, 0.0)
        
        # Craqueo secundario limitado por el tiempo de residencia con gas de arrastre N2
        # Con flujo de N2 (m_dot_N2), los vapores de alquitrán son evacuados rápidamente (tau_res ~ 2-4 s)
        # por lo que solo una fracción menor (~10-25%) del alquitrán alcanza a craquearse a gas en el lecho.
        tau_res_vapor = 3.0 # segundos de tiempo de residencia promedio de volátiles
        frac_cracked = np.clip(k_crk_g * tau_res_vapor, 0.0, 0.35)
        
        r_crk_g = frac_cracked * (self.y_tar_hce * r_hce + self.y_tar_cell * r_cell + self.y_tar_lig * r_lig + self.y_tar_lip * r_lip)
        r_crk_c = 0.08 * r_crk_g # Formación de char secundario por policondensación
        
        # Producción neta de fases
        d_char_dt = (
            self.y_char_hce * r_hce +
            self.y_char_cell * r_cell +
            self.y_char_lig * r_lig +
            self.y_char_lip * r_lip +
            r_crk_c  # Char secundario
        )
        
        # Biovolátiles / Bio-oil (alquitranes y condensables arrastrados por N2 al condensador inferior)
        d_tar_dt = (
            self.y_tar_hce * r_hce +
            self.y_tar_cell * r_cell +
            self.y_tar_lig * r_lig +
            self.y_tar_lip * r_lip -
            r_crk_g - r_crk_c
        )
        d_tar_dt = np.maximum(d_tar_dt, 0.0)
        
        # Biogás (gases permanentes no condensables: CO, CO2, CH4, H2)
        d_gas_dt = (
            self.y_gas_hce * r_hce +
            self.y_gas_cell * r_cell +
            self.y_gas_lig * r_lig +
            self.y_gas_lip * r_lip +
            r_crk_g
        )
        
        d_vap_dt = r_evap
        
        # Término fuente térmico volumétrico neto (W/m^3)
        # q_dot > 0 es exotérmico, q_dot < 0 es endotérmico
        q_dot_rxn = (
            r_evap * self.dH_evap +
            r_hce  * self.dH_hce +
            r_cell * self.dH_cell +
            r_lig  * self.dH_lig +
            r_lip  * self.dH_lip +
            r_crk_g * self.dH_crk_gas +
            r_crk_c * self.dH_crk_char
        )
        
        rates = {
            'r_evap': r_evap,
            'r_hce': r_hce,
            'r_cell': r_cell,
            'r_lig': r_lig,
            'r_lip': r_lip,
            'r_crk_g': r_crk_g,
            'r_crk_c': r_crk_c,
            'd_char_dt': d_char_dt,
            'd_tar_dt': d_tar_dt,
            'd_gas_dt': d_gas_dt,
            'd_vap_dt': d_vap_dt,
            'q_dot_rxn': q_dot_rxn
        }
        
        return rates
