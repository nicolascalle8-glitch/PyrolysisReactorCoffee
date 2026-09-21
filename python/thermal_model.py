"""
thermal_model.py
Modelo bidimensional axisimétrico (r, z) de transferencia de calor y masa en lecho cilíndrico poroso.
Resuelve la ecuación diferencial parcial parabólica no lineal con coeficientes variables:

  rho*Cp * dT/dt = (1/r)*d/dr(r*k*dT/dr) + d/dz(k*dT/dz) - rho_g*Cp_g*v_g*dT/dz + q_dot_rxn

Incluye:
- Discretización espacial en diferencias finitas centradas
- Tratamiento de la singularidad matemática en el centro r = 0 mediante la regla de L'Hôpital:
      lim_{r->0} (1/r)*dT/dr = d^2T/dr^2  =>  Conducción en r=0: 4*k*(T_1 - T_0)/dr^2
- Condiciones de frontera de Robin (convección y radiación) en r = R
- Convección interna por flujo de Darcy de gases pirolíticos (v_g)
- Acoplamiento término a término con la cinética de pirólisis (kinetics.py)
"""

import numpy as np

class ReactorThermal2D:
    def __init__(self, geo_cfg, bio_cfg, therm_cfg, kinetics):
        self.geo = geo_cfg
        self.bio = bio_cfg
        self.therm = therm_cfg
        self.kin = kinetics
        
        self.Nr = self.geo.Nr
        self.Nz = self.geo.Nz
        self.dr = self.geo.dr
        self.dz = self.geo.dz
        self.r = self.geo.r_grid
        self.z = self.geo.z_grid
        
        # Estado inicial de concentraciones másicas (kg/m^3 de lecho)
        self.rho_w0    = self.bio.rho_bulk * self.bio.w_H2O
        self.rho_cell0 = self.bio.rho_bulk * self.bio.w_CELL
        self.rho_hce0  = self.bio.rho_bulk * self.bio.w_HCE
        self.rho_lig0  = self.bio.rho_bulk * self.bio.w_LIG
        self.rho_lip0  = self.bio.rho_bulk * self.bio.w_LIP
        self.rho_char0 = self.bio.rho_bulk * self.bio.w_ASH # Cenizas iniciales
        self.rho_tar0  = 0.0
        self.rho_gas0  = 0.0

    def init_state_arrays(self, T_init=298.15):
        """Inicializa las matrices bidimensionales de temperatura y especies (Nr, Nz)"""
        T = np.full((self.Nr, self.Nz), T_init)
        rho_w    = np.full((self.Nr, self.Nz), self.rho_w0)
        rho_cell = np.full((self.Nr, self.Nz), self.rho_cell0)
        rho_hce  = np.full((self.Nr, self.Nz), self.rho_hce0)
        rho_lig  = np.full((self.Nr, self.Nz), self.rho_lig0)
        rho_lip  = np.full((self.Nr, self.Nz), self.rho_lip0)
        rho_char = np.full((self.Nr, self.Nz), self.rho_char0)
        rho_tar  = np.full((self.Nr, self.Nz), self.rho_tar0)
        rho_gas  = np.full((self.Nr, self.Nz), self.rho_gas0)
        
        return {
            'T': T,
            'w': rho_w,
            'cell': rho_cell,
            'hce': rho_hce,
            'lig': rho_lig,
            'lip': rho_lip,
            'char': rho_char,
            'tar': rho_tar,
            'gas': rho_gas
        }

    def compute_thermophysical_props(self, T, rho_solid):
        """
        Calcula densidad aparente efectiva, capacidad calorífica Cp(T) y conductividad k(T)
        en cada celda de la malla 2D.
        """
        # Fracción de conversión sólida X = 1 - rho_solid / rho_solid0
        rho_solid0 = self.bio.rho_bulk
        conv = np.clip(1.0 - rho_solid / rho_solid0, 0.0, 1.0)
        
        # Conductividad efectiva ponderada (W/(m·K))
        # Aumenta con la temperatura y la formación de biochar poroso
        k_eff = (1.0 - conv) * (self.therm.k_biomass_0 + self.therm.k_biomass_slope * (T - 298.15)) + \
                conv * (self.therm.k_char_0 + self.therm.k_char_slope * (T - 298.15))
        k_eff = np.maximum(k_eff, 0.05)
        
        # Capacidad calorífica específica efectiva Cp (J/(kg·K))
        # Biomasa virgen: ~ 1200 + 1.8*(T - 273)
        # Biochar: ~ 700 + 2.1*(T - 273)
        cp_biomass = 1250.0 + 1.8 * (T - 273.15)
        cp_char = 720.0 + 2.2 * (T - 273.15)
        cp_eff = (1.0 - conv) * cp_biomass + conv * cp_char
        cp_eff = np.maximum(cp_eff, 500.0)
        
        # Densidad volumétrica actual de masa
        rho_eff = np.maximum(rho_solid, 50.0)
        
        return rho_eff, cp_eff, k_eff

    def compute_derivatives(self, state, T_furnace, mode='robin'):
        """
        Calcula las derivadas temporales dT/dt y drho_i/dt en toda la malla 2D (Nr, Nz).
        Implementa rigurosamente la condición en r = 0 derivada en el documento físico.
        """
        T = state['T']
        rho_w    = state['w']
        rho_cell = state['cell']
        rho_hce  = state['hce']
        rho_lig  = state['lig']
        rho_lip  = state['lip']
        rho_char = state['char']
        rho_tar  = state['tar']
        rho_gas  = state['gas']
        
        rho_solid = rho_w + rho_cell + rho_hce + rho_lig + rho_lip + rho_char
        
        # 1. Propiedades térmicas locales
        rho_eff, cp_eff, k_eff = self.compute_thermophysical_props(T, rho_solid)
        
        # 2. Cinética química local en cada nodo
        kin_rates = self.kin.compute_rates(
            T, rho_w, rho_hce, rho_cell, rho_lig, rho_lip, rho_tar, rho_char
        )
        
        # 3. Flujo convectivo de volátiles (flujo de Darcy axial hacia la salida inferior z = L)
        # Estimación de velocidad superficial de gases v_g (m/s) debida a la generación de vapor y gas
        d_volatiles_dt = kin_rates['d_tar_dt'] + kin_rates['d_gas_dt'] + kin_rates['d_vap_dt']
        # Acumulación de caudal hacia abajo: v_g(z) ~ integral_0^z (d_vol/dt)/rho_g dz
        v_g = np.zeros_like(T)
        for j in range(1, self.Nz):
            v_g[:, j] = v_g[:, j-1] + (d_volatiles_dt[:, j] / self.bio.rho_gas) * self.dz
        v_g = np.clip(v_g, 0.0, 0.5) # Limitar a rangos físicos en lechos empacados
        
        # 4. Derivadas de conducción térmica espacial
        dT_dt = np.zeros_like(T)
        dr2 = self.dr ** 2
        dz2 = self.dz ** 2
        
        # ----------------------------------------------------------------------
        # A. NODOS INTERNOS (1 <= i <= Nr-2, 1 <= j <= Nz-2)
        # ----------------------------------------------------------------------
        # Conducción radial: (1/r)*d/dr(r*k*dT/dr) = k*d2T/dr2 + (k/r + dk/dr)*dT/dr
        for i in range(1, self.Nr - 1):
            r_i = self.r[i]
            for j in range(1, self.Nz - 1):
                # Derivada segunda radial
                d2T_dr2 = (T[i+1, j] - 2.0 * T[i, j] + T[i-1, j]) / dr2
                # Derivada primera radial (centrada)
                dT_dr = (T[i+1, j] - T[i-1, j]) / (2.0 * self.dr)
                # Gradiente de conductividad dk/dr
                dk_dr = (k_eff[i+1, j] - k_eff[i-1, j]) / (2.0 * self.dr)
                
                term_radial = k_eff[i, j] * (d2T_dr2 + (1.0 / r_i) * dT_dr) + dk_dr * dT_dr
                
                # Conducción axial
                d2T_dz2 = (T[i, j+1] - 2.0 * T[i, j] + T[i, j-1]) / dz2
                term_axial = k_eff[i, j] * d2T_dz2
                
                # Convección axial de gases pirolíticos (upwind/retrograda en dirección del flujo)
                dT_dz = (T[i, j] - T[i, j-1]) / self.dz
                term_conv = self.bio.rho_gas * self.therm.cp_steam * v_g[i, j] * dT_dz
                
                # Fuente por reacción química
                q_rxn = kin_rates['q_dot_rxn'][i, j]
                
                # Balance neto de energía
                dT_dt[i, j] = (term_radial + term_axial - term_conv + q_rxn) / (rho_eff[i, j] * cp_eff[i, j])

        # ----------------------------------------------------------------------
        # B. CONDICIÓN EN EL CENTRO (r = 0, i = 0): REGLA DE L'HÔPITAL
        # Según la deducción en el documento (Página 3):
        # lim_{r->0} (1/r)*dT/dr = d2T/dr2 => d/dr(r*k*dT/dr)/r = 2*k*d2T/dr2
        # Por simetría dT/dr = 0 => T_{-1} = T_1 => d2T/dr2 = 2*(T_1 - T_0)/dr^2
        # Por tanto: término radial = 4 * k * (T_1 - T_0) / dr^2
        # ----------------------------------------------------------------------
        for j in range(1, self.Nz - 1):
            term_radial_center = 4.0 * k_eff[0, j] * (T[1, j] - T[0, j]) / dr2
            d2T_dz2 = (T[0, j+1] - 2.0 * T[0, j] + T[0, j-1]) / dz2
            term_axial = k_eff[0, j] * d2T_dz2
            
            dT_dz = (T[0, j] - T[0, j-1]) / self.dz
            term_conv = self.bio.rho_gas * self.therm.cp_steam * v_g[0, j] * dT_dz
            q_rxn = kin_rates['q_dot_rxn'][0, j]
            
            dT_dt[0, j] = (term_radial_center + term_axial - term_conv + q_rxn) / (rho_eff[0, j] * cp_eff[0, j])

        # ----------------------------------------------------------------------
        # C. CONDICIÓN EN LA PARED EXTERIOR (r = R, i = Nr - 1)
        # ----------------------------------------------------------------------
        if mode == 'dirichlet':
            # La pared del reactor sigue exactamente la temperatura del horno/rampa
            dT_dt[self.Nr - 1, :] = 0.0
            T[self.Nr - 1, :] = T_furnace
        else:
            # Condición de Robin: flujo de calor por convección y radiación desde el horno
            # -k * dT/dr = h_comb * (T_wall - T_furnace)
            # En diferencias finitas con nodo fantasma: T_{Nr} = T_{Nr-2} + 2*dr*(h/k)*(T_furnace - T_wall)
            h_comb = self.therm.h_furnace_tube
            for j in range(1, self.Nz - 1):
                T_wall = T[self.Nr - 1, j]
                k_w = k_eff[self.Nr - 1, j]
                q_in = h_comb * (T_furnace - T_wall)
                
                # Derivada segunda considerando el flujo entrante
                # dT/dr en la pared = q_in / k_w
                d2T_dr2 = 2.0 * (T[self.Nr - 2, j] - T_wall + (self.dr * q_in / k_w)) / dr2
                r_wall = self.r[-1]
                term_radial = k_w * (d2T_dr2 + (1.0 / r_wall) * (q_in / k_w))
                
                d2T_dz2 = (T[self.Nr - 1, j+1] - 2.0 * T_wall + T[self.Nr - 1, j-1]) / dz2
                term_axial = k_w * d2T_dz2
                q_rxn = kin_rates['q_dot_rxn'][self.Nr - 1, j]
                
                dT_dt[self.Nr - 1, j] = (term_radial + term_axial + q_rxn) / (rho_eff[self.Nr - 1, j] * cp_eff[self.Nr - 1, j])

        # ----------------------------------------------------------------------
        # D. EXTREMOS AXIALES: z = 0 (Entrada N2) y z = L (Salida de vapores)
        # ----------------------------------------------------------------------
        # z = 0: convección con gas N2 entrante a T_amb
        dT_dt[:, 0] = dT_dt[:, 1]
        # z = L: flujo libre / convección de salida
        dT_dt[:, self.Nz - 1] = dT_dt[:, self.Nz - 2]
        
        # 5. Derivadas temporales de concentración de componentes
        d_rho = {
            'w':    -kin_rates['r_evap'],
            'cell': -kin_rates['r_cell'],
            'hce':  -kin_rates['r_hce'],
            'lig':  -kin_rates['r_lig'],
            'lip':  -kin_rates['r_lip'],
            'char':  kin_rates['d_char_dt'],
            'tar':   kin_rates['d_tar_dt'],
            'gas':   kin_rates['d_gas_dt']
        }
        
        return dT_dt, d_rho, kin_rates
