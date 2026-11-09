r"""
kinetics_lumped.py
Modelo Cinético Agrupado (Lumped Model) para la pirólisis de residuos de café (SCG).

Esquema de reacciones en paralelo y consecutivas (Modelo de Shafizadeh-Chin / Thurner & Mann / Di Blasi):

          k1 (Char)     ---> Biochar primario
         /
Biomasa --- k2 (Vol)      ---> Biovolátiles / Bio-oil (Vapores condensables)
         \                                  |
          k3 (GNC)      ---> Gases No Condensables (GNC)
                                            |
         +----------------------------------+
         |
         |---> k4 (Repolimerización) ---> Biochar secundario + GNC
         |---> k5 (Craqueo térmico)  ---> Gases No Condensables (GNC adicionales)

Humedad (H2O_l) --- k_evap ---> Vapor de agua (H2O_v)
"""

import numpy as np

R_UNIVERSAL = 8.31446  # J/(mol·K)

class LumpedPyrolysisKinetics:
    def __init__(self):
        # ----------------------------------------------------------------------
        # PARÁMETROS CINÉTICOS ARRHENIUS: k = A * exp(-Ea / (R*T))
        # Basados en el modelo de tres vías en paralelo de Chan / Di Blasi adaptado a café:
        # Permite reproducir la partición experimental (~30-35% Biochar, ~40-45% Volátiles, ~18-22% GNC)
        # ----------------------------------------------------------------------
        # Reacción 1: Biomasa -> Biochar primario
        self.A1 = 1.15e8        # s^-1
        self.Ea1 = 122.0e3      # J/mol
        self.dH1 = 120.0e3      # J/kg (exotérmico)

        # Reacción 2: Biomasa -> Biovolátiles / Bio-oil (Vapores condensables)
        self.A2 = 1.65e8        # s^-1
        self.Ea2 = 122.5e3      # J/mol
        self.dH2 = -210.0e3     # J/kg (endotérmico)

        # Reacción 3: Biomasa -> Gases No Condensables primarios (GNC: CO, CO2, CH4, H2)
        self.A3 = 7.50e7        # s^-1
        self.Ea3 = 118.0e3      # J/mol
        self.dH3 = -80.0e3      # J/kg (endotérmico)

        # Reacción 4: Biovolátiles -> Biochar secundario (repolimerización heterogénea)
        self.A4 = 1.0e5         # s^-1
        self.Ea4 = 95.0e3       # J/mol
        self.dH4 = 60.0e3       # J/kg (exotérmico)

        # Reacción 5: Biovolátiles -> Gases No Condensables (craqueo térmico homogéneo)
        self.A5 = 4.28e6        # s^-1
        self.Ea5 = 108.0e3      # J/mol
        self.dH5 = -42.0e3      # J/kg

        # Evaporación de humedad
        self.A_evap = 5.0e5     # s^-1
        self.Ea_evap = 41.8e3   # J/mol
        self.dH_evap = -2.26e6  # J/kg (fuertemente endotérmico)

        # Tiempo de residencia promedio de vapores en el lecho con purga de N2 (s)
        self.tau_res_vapor = 3.0

    def compute_rates(self, T, rho_biomass, rho_water, rho_volatiles, rho_char):
        """
        Calcula las velocidades de consumo y producción de cada pseudocomponente agrupado (LUMPED).
        
        Parámetros:
          T: Temperatura local (K)
          rho_biomass: Concentración de biomasa virgen de café (kg/m^3)
          rho_water: Concentración de humedad líquida (kg/m^3)
          rho_volatiles: Concentración de biovolátiles/alquitranes presentes (kg/m^3)
          rho_char: Concentración de biochar presente (kg/m^3)
        """
        T_clipped = np.maximum(T, 273.15)
        RT = R_UNIVERSAL * T_clipped

        # 1. Constantes de velocidad de Arrhenius (1/s)
        k1 = self.A1 * np.exp(-self.Ea1 / RT)
        k2 = self.A2 * np.exp(-self.Ea2 / RT)
        k3 = self.A3 * np.exp(-self.Ea3 / RT)
        k4 = self.A4 * np.exp(-self.Ea4 / RT)
        k5 = self.A5 * np.exp(-self.Ea5 / RT)

        k_evap = np.where(T_clipped > 350.0, self.A_evap * np.exp(-self.Ea_evap / RT), 0.0)

        # 2. Tasas primarias de descomposición de la biomasa (kg/(m^3·s))
        rho_B_pos = np.maximum(rho_biomass, 0.0)
        r1_char_prim = k1 * rho_B_pos       # Biomasa -> Biochar
        r2_vol_prim  = k2 * rho_B_pos       # Biomasa -> Biovolátiles
        r3_gnc_prim  = k3 * rho_B_pos       # Biomasa -> Gases No Condensables (GNC)

        # Tasa total de consumo de biomasa de café
        r_biomass_decay = r1_char_prim + r2_vol_prim + r3_gnc_prim

        # Evaporación de agua
        r_evap = k_evap * np.maximum(rho_water, 0.0)

        # 3. Reacciones secundarias de los biovolátiles (Tar Cracking y Repolimerización)
        # Modulado por el arrastre de gas N2: solo una fracción permanece el tiempo suficiente para craquear
        frac_craqueo = np.clip(k5 * self.tau_res_vapor, 0.0, 0.30)
        r5_gnc_sec  = frac_craqueo * r2_vol_prim
        r4_char_sec = 0.08 * r5_gnc_sec     # Formación de char secundario

        # 4. Derivadas netas de producción de cada fase (kg/(m^3·s))
        # A) Consumo de Biomasa
        d_biomass_dt = -r_biomass_decay

        # B) Producción neta de Biochar (primario + secundario)
        d_char_dt = r1_char_prim + r4_char_sec

        # C) Producción neta de Biovolátiles / Bio-oil (que salen al condensador)
        d_volatiles_dt = np.maximum(r2_vol_prim - r4_char_sec - r5_gnc_sec, 0.0)

        # D) Producción neta de Gases No Condensables (GNC: primarios + secundarios)
        d_gnc_dt = r3_gnc_prim + r5_gnc_sec

        # E) Vapor de agua
        d_vap_dt = r_evap

        # 5. Término fuente de calor volumétrico q_rxn (W/m^3)
        q_dot_rxn = (
            r_evap       * self.dH_evap +
            r1_char_prim * self.dH1 +
            r2_vol_prim  * self.dH2 +
            r3_gnc_prim  * self.dH3 +
            r4_char_sec  * self.dH4 +
            r5_gnc_sec   * self.dH5
        )

        return {
            'd_biomass_dt': d_biomass_dt,
            'd_char_dt': d_char_dt,
            'd_volatiles_dt': d_volatiles_dt,
            'd_gnc_dt': d_gnc_dt,
            'd_vap_dt': d_vap_dt,
            'r_evap': r_evap,
            'r_biomass_decay': r_biomass_decay,
            'r_char_prim': r1_char_prim,
            'r_vol_prim': r2_vol_prim,
            'r_gnc_prim': r3_gnc_prim,
            'r_char_sec': r4_char_sec,
            'r_gnc_sec': r5_gnc_sec,
            'q_dot_rxn': q_dot_rxn
        }
