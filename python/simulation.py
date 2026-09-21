"""
simulation.py
Motor de integración temporal para el reactor vertical de pirólisis.
Acopla la cinética, la transferencia de calor 2D y el control PID del horno eléctrico.
Calcula la evolución del mapa de calor, masa residual y rendimientos porcentuales
de Biochar, Biogás y Biovolátiles (Bio-oil/Tar).
"""

import numpy as np
from config import GeometryConfig, CoffeeBiomassConfig, ThermalPropertiesConfig, FurnaceConfig, HeatingRampConfig
from kinetics import PyrolysisKinetics
from thermal_model import ReactorThermal2D
from furnace_control import ElectricFurnaceModel

class PyrolysisSimulator:
    def __init__(
        self,
        geo_cfg=None,
        bio_cfg=None,
        therm_cfg=None,
        furnace_cfg=None,
        ramp_cfg=None
    ):
        self.geo = geo_cfg if geo_cfg is not None else GeometryConfig()
        self.bio = bio_cfg if bio_cfg is not None else CoffeeBiomassConfig()
        self.therm = therm_cfg if therm_cfg is not None else ThermalPropertiesConfig()
        self.furnace_cfg = furnace_cfg if furnace_cfg is not None else FurnaceConfig()
        self.ramp = ramp_cfg if ramp_cfg is not None else HeatingRampConfig()
        
        self.kin = PyrolysisKinetics()
        self.reactor = ReactorThermal2D(self.geo, self.bio, self.therm, self.kin)
        self.furnace = ElectricFurnaceModel(self.furnace_cfg, self.ramp)
        
        # Volúmenes de control de cada anillo cilíndrico (2*pi*r*dr*dz)
        self.dV = np.zeros((self.geo.Nr, self.geo.Nz))
        dr = self.geo.dr
        dz = self.geo.dz
        for i in range(self.geo.Nr):
            r_i = self.geo.r_grid[i]
            if i == 0:
                # Cilindro central (r de 0 a dr/2)
                area_annulus = np.pi * (dr / 2.0)**2
            elif i == self.geo.Nr - 1:
                # Anillo exterior
                area_annulus = np.pi * (self.geo.R**2 - (self.geo.R - dr/2.0)**2)
            else:
                area_annulus = 2.0 * np.pi * r_i * dr
            self.dV[i, :] = area_annulus * dz
            
        self.V_total = np.sum(self.dV)

    def integrate_mass(self, rho_matrix):
        """Calcula la masa total integrada (kg) de un componente en el reactor"""
        return np.sum(rho_matrix * self.dV)

    def run(self, max_dt=0.5, save_interval_s=15.0, verbose=True):
        """
        Ejecuta la simulación numérica completa a lo largo del tiempo de rampa y meseta.
        
        max_dt: paso temporal máximo (s) que respeta el criterio CFL de difusión térmica
        save_interval_s: intervalo de tiempo para registrar datos de visualización
        """
        t_final = self.ramp.total_time
        state = self.reactor.init_state_arrays(T_init=self.ramp.T_init_K)
        
        # Masa inicial total de biomasa empacada
        M_init_total = (
            self.integrate_mass(state['w']) +
            self.integrate_mass(state['cell']) +
            self.integrate_mass(state['hce']) +
            self.integrate_mass(state['lig']) +
            self.integrate_mass(state['lip']) +
            self.integrate_mass(state['char'])
        )
        
        if verbose:
            print(f"===========================================================")
            print(f"  SIMULADOR DE PIRÓLISIS DE RESIDUOS DE CAFÉ (SCG)")
            print(f"===========================================================")
            print(f"Geometría: R = {self.geo.R*1000:.1f} mm, L = {self.geo.L*100:.1f} cm")
            print(f"Masa inicial de café empacado: {M_init_total*1000:.2f} g")
            print(f"Humedad inicial: {self.bio.w_H2O*100:.1f} % wt")
            print(f"Rampa: {self.ramp.ramp_rate_K_s*60:.1f} °C/min hasta {(self.ramp.T_final_K-273.15):.1f} °C")
            print(f"Tiempo total: {t_final/60.0:.1f} min ({t_final:.0f} s)")
            print(f"Malla: Nr = {self.geo.Nr}, Nz = {self.geo.Nz}")
            print(f"-----------------------------------------------------------")

        # Acumuladores de masa de productos gaseosos que abandonan el lecho
        cum_gas_produced = 0.0
        cum_tar_produced = 0.0
        cum_vap_produced = 0.0
        
        # Historial de registros
        time_history = []
        T_center_history = []
        T_mid_history = []
        T_wall_history = []
        T_furnace_history = []
        T_setpoint_history = []
        residual_mass_fraction = []
        
        # Rendimientos acumulados (%)
        yield_biochar_history = []
        yield_biogas_history = []
        yield_biovolatiles_history = []
        yield_water_history = []
        
        # Mapas de calor 2D guardados en instantes clave
        snapshots = []
        snapshot_times = [
            0.0,
            0.30 * self.ramp.t_ramp_end,
            0.60 * self.ramp.t_ramp_end,
            self.ramp.t_ramp_end,
            t_final
        ]
        next_snap_idx = 0

        t = 0.0
        last_save_t = -save_interval_s
        mid_r_idx = self.geo.Nr // 2
        mid_z_idx = self.geo.Nz // 2

        while t <= t_final:
            # 1. Control del horno eléctrico y temperatura de pared
            T_wall_avg = np.mean(state['T'][-1, :])
            T_furnace, T_set, duty = self.furnace.step(t, max_dt, T_wall_avg)
            
            # 2. Paso de CFL adaptativo para estabilidad
            # C_diff = alpha * dt / dr^2 <= 0.4
            alpha_approx = self.therm.k_biomass_0 / (self.bio.rho_bulk * 1500.0)
            cfl_dt = 0.35 * (self.geo.dr ** 2) / max(alpha_approx, 1e-8)
            dt_step = min(max_dt, cfl_dt, t_final - t if t_final > t else max_dt)
            if dt_step <= 0:
                break
                
            # 3. Derivadas del reactor
            dT_dt, d_rho, kin_rates = self.reactor.compute_derivatives(state, T_furnace, mode='robin')
            
            # 4. Integración temporal (Heun / predictor-corrector para alta estabilidad)
            # Predictor (Euler)
            T_pred = state['T'] + dT_dt * dt_step
            state_pred = {k: np.maximum(state[k] + d_rho[k] * dt_step, 0.0) for k in d_rho}
            state_pred['T'] = T_pred
            
            # Corrector
            dT_dt_corr, d_rho_corr, kin_rates_corr = self.reactor.compute_derivatives(state_pred, T_furnace, mode='robin')
            
            state['T'] += 0.5 * (dT_dt + dT_dt_corr) * dt_step
            for k in d_rho:
                state[k] = np.maximum(state[k] + 0.5 * (d_rho[k] + d_rho_corr[k]) * dt_step, 0.0)
                
            # 5. Acumular gases y vapores liberados
            dV_matrix = self.dV
            d_gas_total = np.sum(0.5 * (kin_rates['d_gas_dt'] + kin_rates_corr['d_gas_dt']) * dV_matrix) * dt_step
            d_tar_total = np.sum(0.5 * (kin_rates['d_tar_dt'] + kin_rates_corr['d_tar_dt']) * dV_matrix) * dt_step
            d_vap_total = np.sum(0.5 * (kin_rates['d_vap_dt'] + kin_rates_corr['d_vap_dt']) * dV_matrix) * dt_step
            
            cum_gas_produced += d_gas_total
            cum_tar_produced += d_tar_total
            cum_vap_produced += d_vap_total
            
            # 6. Registrar snapshots para el mapa de calor
            if next_snap_idx < len(snapshot_times) and t >= snapshot_times[next_snap_idx]:
                snapshots.append({
                    'time': t,
                    'T_map': np.copy(state['T']),
                    'char_map': np.copy(state['char']),
                    'label': f"t = {t/60.0:.1f} min"
                })
                next_snap_idx += 1

            # 7. Guardar historial periódico
            if (t - last_save_t) >= save_interval_s or t >= t_final:
                last_save_t = t
                time_history.append(t)
                T_center_history.append(state['T'][0, mid_z_idx])
                T_mid_history.append(state['T'][mid_r_idx, mid_z_idx])
                T_wall_history.append(state['T'][-1, mid_z_idx])
                T_furnace_history.append(T_furnace)
                T_setpoint_history.append(T_set)
                
                # Masas presentes
                M_char = self.integrate_mass(state['char'])
                M_solid_remaining = (
                    self.integrate_mass(state['w']) +
                    self.integrate_mass(state['cell']) +
                    self.integrate_mass(state['hce']) +
                    self.integrate_mass(state['lig']) +
                    self.integrate_mass(state['lip']) +
                    M_char
                )
                residual_mass_fraction.append((M_solid_remaining / M_init_total) * 100.0)
                
                yield_biochar_history.append((M_char / M_init_total) * 100.0)
                yield_biogas_history.append((cum_gas_produced / M_init_total) * 100.0)
                yield_biovolatiles_history.append((cum_tar_produced / M_init_total) * 100.0)
                yield_water_history.append((cum_vap_produced / M_init_total) * 100.0)
                
                if verbose and (len(time_history) % 15 == 0 or t >= t_final):
                    print(
                        f"t: {t/60.0:5.1f} min | "
                        f"T_horno: {T_furnace-273.15:5.1f}°C | "
                        f"T_pared: {state['T'][-1, mid_z_idx]-273.15:5.1f}°C | "
                        f"T_centro: {state['T'][0, mid_z_idx]-273.15:5.1f}°C | "
                        f"Biochar: {yield_biochar_history[-1]:5.1f}% | "
                        f"Biogás: {yield_biogas_history[-1]:5.1f}% | "
                        f"Biovolátiles: {yield_biovolatiles_history[-1]:5.1f}%"
                    )

            t += dt_step

        # Asegurar último snapshot
        if len(snapshots) < len(snapshot_times):
            snapshots.append({
                'time': t_final,
                'T_map': np.copy(state['T']),
                'char_map': np.copy(state['char']),
                'label': f"Final (t = {t_final/60.0:.1f} min)"
            })

        M_unconverted = (
            self.integrate_mass(state['cell']) +
            self.integrate_mass(state['hce']) +
            self.integrate_mass(state['lig']) +
            self.integrate_mass(state['lip'])
        )
        unconverted_pct = (M_unconverted / M_init_total) * 100.0

        final_yields = {
            'biochar_pct': yield_biochar_history[-1],
            'biogas_pct': yield_biogas_history[-1],
            'biovolatiles_pct': yield_biovolatiles_history[-1],
            'water_vapor_pct': yield_water_history[-1],
            'unconverted_pct': unconverted_pct
        }
        
        balance_check = sum(final_yields.values())
        if verbose:
            print(f"-----------------------------------------------------------")
            print(f"RENDIMIENTOS FINALES:")
            print(f"  • Biochar (sólido carbonoso): {final_yields['biochar_pct']:.2f} %")
            print(f"  • Biogás (gases no condensables): {final_yields['biogas_pct']:.2f} %")
            print(f"  • Biovolátiles / Bio-oil (alquitranes/lípidos): {final_yields['biovolatiles_pct']:.2f} %")
            print(f"  • Vapor de agua (humedad + agua de pirolisis): {final_yields['water_vapor_pct']:.2f} %")
            if unconverted_pct > 0.1:
                print(f"  • Biomasa no convertida: {unconverted_pct:.2f} %")
            print(f"  TOTAL BALANCE DE MASA: {balance_check:.2f} %")
            print(f"===========================================================")

        results = {
            'time_s': np.array(time_history),
            'time_min': np.array(time_history) / 60.0,
            'T_center': np.array(T_center_history) - 273.15,
            'T_mid': np.array(T_mid_history) - 273.15,
            'T_wall': np.array(T_wall_history) - 273.15,
            'T_furnace': np.array(T_furnace_history) - 273.15,
            'T_setpoint': np.array(T_setpoint_history) - 273.15,
            'residual_mass_pct': np.array(residual_mass_fraction),
            'yield_biochar': np.array(yield_biochar_history),
            'yield_biogas': np.array(yield_biogas_history),
            'yield_biovolatiles': np.array(yield_biovolatiles_history),
            'yield_water': np.array(yield_water_history),
            'final_yields': final_yields,
            'snapshots': snapshots,
            'geo': self.geo,
            'ramp': self.ramp
        }
        
        return results
