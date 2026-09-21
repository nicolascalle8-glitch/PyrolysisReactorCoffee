"""
thermal_model_lumped.py
Modelo bidimensional axisimétrico (r, z) acoplado al modelo cinético agrupado (Lumped).
Especies de estado:
1. T: Temperatura local (K)
2. biomass: Concentración de biomasa virgen de café (kg/m^3)
3. char: Concentración de biochar (kg/m^3)
4. water: Concentración de humedad líquida (kg/m^3)
"""

import numpy as np

class ReactorThermal2DLumped:
    def __init__(self, geo_cfg, bio_cfg, therm_cfg, kinetics_lumped):
        self.geo = geo_cfg
        self.bio = bio_cfg
        self.therm = therm_cfg
        self.kin = kinetics_lumped
        
        self.Nr = self.geo.Nr
        self.Nz = self.geo.Nz
        self.dr = self.geo.dr
        self.dz = self.geo.dz
        self.r = self.geo.r_grid
        self.z = self.geo.z_grid
        
        # Concentraciones iniciales
        self.rho_water_0   = self.bio.rho_bulk * self.bio.w_H2O
        self.rho_char_0    = self.bio.rho_bulk * self.bio.w_ASH # Cenizas iniciales
        # Biomasa orgánica seca inicial
        self.rho_biomass_0 = self.bio.rho_bulk - self.rho_water_0 - self.rho_char_0

    def init_state(self, T_init=298.15):
        """Inicializa las matrices de estado 2D (Nr, Nz)"""
        T = np.full((self.Nr, self.Nz), T_init)
        biomass = np.full((self.Nr, self.Nz), self.rho_biomass_0)
        char    = np.full((self.Nr, self.Nz), self.rho_char_0)
        water   = np.full((self.Nr, self.Nz), self.rho_water_0)
        
        return {
            'T': T,
            'biomass': biomass,
            'char': char,
            'water': water
        }

    def compute_derivatives(self, state, T_furnace):
        """
        Calcula las derivadas temporales dT/dt, d(biomass)/dt, d(char)/dt y d(water)/dt.
        Aplica la regla de L'Hôpital en r = 0 y la condición de Robin en r = R.
        """
        T = state['T']
        biomass = state['biomass']
        char = state['char']
        water = state['water']
        
        rho_solid = biomass + char + water
        
        # 1. Grado de conversión X
        conv = np.clip(1.0 - (biomass + water) / (self.rho_biomass_0 + self.rho_water_0), 0.0, 1.0)
        
        # Conductividad y Cp variables
        k_eff = (1.0 - conv) * (self.therm.k_biomass_0 + self.therm.k_biomass_slope * (T - 298.15)) + \
                conv * (self.therm.k_char_0 + self.therm.k_char_slope * (T - 298.15))
        k_eff = np.maximum(k_eff, 0.05)
        
        cp_eff = (1.0 - conv) * (1250.0 + 1.8 * (T - 273.15)) + conv * (720.0 + 2.2 * (T - 273.15))
        cp_eff = np.maximum(cp_eff, 500.0)
        
        rho_eff = np.maximum(rho_solid, 50.0)
        
        # 2. Cinética Lumped
        rates = self.kin.compute_rates(
            T,
            rho_biomass=biomass,
            rho_water=water,
            rho_volatiles=np.zeros_like(biomass),
            rho_char=char
        )
        
        # 3. Flujo convectivo de volátiles (Darcy) hacia la salida inferior z = L
        d_vol_dt = rates['d_volatiles_dt'] + rates['d_gnc_dt'] + rates['d_vap_dt']
        v_g = np.zeros_like(T)
        for j in range(1, self.Nz):
            v_g[:, j] = v_g[:, j-1] + (d_vol_dt[:, j] / self.bio.rho_gas) * self.dz
        v_g = np.clip(v_g, 0.0, 0.5)
        
        # 4. Derivadas de conducción de calor
        dT_dt = np.zeros_like(T)
        dr2 = self.dr ** 2
        dz2 = self.dz ** 2
        
        # A. Nodos internos (1 <= i <= Nr-2, 1 <= j <= Nz-2)
        for i in range(1, self.Nr - 1):
            r_i = self.r[i]
            for j in range(1, self.Nz - 1):
                d2T_dr2 = (T[i+1, j] - 2.0 * T[i, j] + T[i-1, j]) / dr2
                dT_dr   = (T[i+1, j] - T[i-1, j]) / (2.0 * self.dr)
                dk_dr   = (k_eff[i+1, j] - k_eff[i-1, j]) / (2.0 * self.dr)
                term_radial = k_eff[i, j] * (d2T_dr2 + (1.0 / r_i) * dT_dr) + dk_dr * dT_dr
                
                d2T_dz2 = (T[i, j+1] - 2.0 * T[i, j] + T[i, j-1]) / dz2
                term_axial = k_eff[i, j] * d2T_dz2
                
                dT_dz = (T[i, j] - T[i, j-1]) / self.dz
                term_conv = self.bio.rho_gas * self.therm.cp_steam * v_g[i, j] * dT_dz
                
                q_rxn = rates['q_dot_rxn'][i, j]
                dT_dt[i, j] = (term_radial + term_axial - term_conv + q_rxn) / (rho_eff[i, j] * cp_eff[i, j])

        # B. Centro r = 0 (Regla de L'Hôpital): 4 * k * (T_1 - T_0) / dr^2
        for j in range(1, self.Nz - 1):
            term_radial_center = 4.0 * k_eff[0, j] * (T[1, j] - T[0, j]) / dr2
            d2T_dz2 = (T[0, j+1] - 2.0 * T[0, j] + T[0, j-1]) / dz2
            term_axial = k_eff[0, j] * d2T_dz2
            
            dT_dz = (T[0, j] - T[0, j-1]) / self.dz
            term_conv = self.bio.rho_gas * self.therm.cp_steam * v_g[0, j] * dT_dz
            q_rxn = rates['q_dot_rxn'][0, j]
            
            dT_dt[0, j] = (term_radial_center + term_axial - term_conv + q_rxn) / (rho_eff[0, j] * cp_eff[0, j])

        # C. Pared r = R (Robin con horno e inercia de piedra)
        h_comb = self.therm.h_furnace_tube
        for j in range(1, self.Nz - 1):
            T_wall = T[self.Nr - 1, j]
            k_w = k_eff[self.Nr - 1, j]
            q_in = h_comb * (T_furnace - T_wall)
            
            d2T_dr2 = 2.0 * (T[self.Nr - 2, j] - T_wall + (self.dr * q_in / k_w)) / dr2
            term_radial = k_w * (d2T_dr2 + (1.0 / self.r[-1]) * (q_in / k_w))
            
            d2T_dz2 = (T[self.Nr - 1, j+1] - 2.0 * T_wall + T[self.Nr - 1, j-1]) / dz2
            term_axial = k_w * d2T_dz2
            q_rxn = rates['q_dot_rxn'][self.Nr - 1, j]
            
            dT_dt[self.Nr - 1, j] = (term_radial + term_axial + q_rxn) / (rho_eff[self.Nr - 1, j] * cp_eff[self.Nr - 1, j])

        # D. Extremos axiales (z = 0 y z = L)
        dT_dt[:, 0] = dT_dt[:, 1]
        dT_dt[:, self.Nz - 1] = dT_dt[:, self.Nz - 2]
        
        # 5. Derivadas de componentes agrupados
        d_biomass = rates['d_biomass_dt']
        d_char    = rates['d_char_dt']
        d_water   = -rates['r_evap']
        
        return dT_dt, d_biomass, d_char, d_water, rates
