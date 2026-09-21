"""
simulation_lumped.py
Motor de integración temporal para el modelo agrupado (Lumped Model).
Calcula la evolución de la temperatura, degradación de biomasa y los rendimientos
de:
1. Biochar (sólido carbonoso)
2. Biovolátiles / Bio-oil (líquidos condensables)
3. Gases No Condensables (GNC: CO, CO2, CH4, H2)
4. Vapor de agua (humedad inicial)
"""

import numpy as np
from config import GeometryConfig, CoffeeBiomassConfig, ThermalPropertiesConfig, FurnaceConfig, HeatingRampConfig
from kinetics_lumped import LumpedPyrolysisKinetics
from thermal_model_lumped import ReactorThermal2DLumped
from furnace_control import ElectricFurnaceModel

class LumpedPyrolysisSimulator:
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
        
        self.kin = LumpedPyrolysisKinetics()
        self.reactor = ReactorThermal2DLumped(self.geo, self.bio, self.therm, self.kin)
        self.furnace = ElectricFurnaceModel(self.furnace_cfg, self.ramp)
        
        # Volúmenes de control diferenciales para coordenadas cilíndricas dV = 2*pi*r*dr*dz
        self.dV = np.zeros((self.geo.Nr, self.geo.Nz))
        dr = self.geo.dr
        dz = self.geo.dz
        for i in range(self.geo.Nr):
            r_i = self.geo.r_grid[i]
            if i == 0:
                area_ring = np.pi * (dr / 2.0)**2
            elif i == self.geo.Nr - 1:
                area_ring = np.pi * (self.geo.R**2 - (self.geo.R - dr/2.0)**2)
            else:
                area_ring = 2.0 * np.pi * r_i * dr
            self.dV[i, :] = area_ring * dz

    def integrate_mass(self, matrix):
        return np.sum(matrix * self.dV)

    def run(self, max_dt=0.5, save_interval_s=15.0, verbose=True):
        t_final = self.ramp.total_time
        state = self.reactor.init_state(T_init=self.ramp.T_init_K)
        
        M_init_total = (
            self.integrate_mass(state['biomass']) +
            self.integrate_mass(state['char']) +
            self.integrate_mass(state['water'])
        )
        
        if verbose:
            print(f"===========================================================")
            print(f"  SIMULADOR DE PIRÓLISIS (MODELO CINÉTICO AGRUPADO - LUMPED)")
            print(f"===========================================================")
            print(f"Masa inicial total: {M_init_total*1000:.2f} g")
            print(f"Humedad: {self.bio.w_H2O*100:.1f} % | Cenizas: {self.bio.w_ASH*100:.2f} %")
            print(f"Rampa: {self.ramp.ramp_rate_K_s*60:.1f} °C/min hasta {(self.ramp.T_final_K-273.15):.1f} °C")
            print(f"Tiempo total: {t_final/60.0:.1f} min ({t_final:.0f} s)")
            print(f"-----------------------------------------------------------")

        cum_volatiles_produced = 0.0
        cum_gnc_produced = 0.0
        cum_water_evaporated = 0.0
        
        time_history = []
        T_center_history = []
        T_mid_history = []
        T_wall_history = []
        T_furnace_history = []
        T_setpoint_history = []
        
        # Historial de porcentajes
        yield_biochar_history = []
        yield_volatiles_history = []
        yield_gnc_history = []
        yield_water_history = []
        residual_mass_fraction = []
        
        # Tasas de reacción promedio en el lecho (para ver claramente el consumo)
        rate_biomass_decay_history = []
        rate_volatiles_history = []
        rate_char_history = []
        rate_gnc_history = []

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
            T_wall_avg = np.mean(state['T'][-1, :])
            T_furnace, T_set, duty = self.furnace.step(t, max_dt, T_wall_avg)
            
            # Criterio CFL difusivo
            alpha_approx = self.therm.k_biomass_0 / (self.bio.rho_bulk * 1500.0)
            cfl_dt = 0.35 * (self.geo.dr ** 2) / max(alpha_approx, 1e-8)
            dt_step = min(max_dt, cfl_dt, t_final - t if t_final > t else max_dt)
            if dt_step <= 0:
                break
                
            # Derivadas (Predictor-Corrector Heun)
            dT_dt, d_bio, d_char, d_water, rates = self.reactor.compute_derivatives(state, T_furnace)
            
            # Predictor
            T_pred = state['T'] + dT_dt * dt_step
            bio_pred = np.maximum(state['biomass'] + d_bio * dt_step, 0.0)
            char_pred = np.maximum(state['char'] + d_char * dt_step, 0.0)
            water_pred = np.maximum(state['water'] + d_water * dt_step, 0.0)
            
            state_pred = {'T': T_pred, 'biomass': bio_pred, 'char': char_pred, 'water': water_pred}
            
            # Corrector
            dT_dt_c, d_bio_c, d_char_c, d_water_c, rates_c = self.reactor.compute_derivatives(state_pred, T_furnace)
            
            state['T'] += 0.5 * (dT_dt + dT_dt_c) * dt_step
            state['biomass'] = np.maximum(state['biomass'] + 0.5 * (d_bio + d_bio_c) * dt_step, 0.0)
            state['char'] = np.maximum(state['char'] + 0.5 * (d_char + d_char_c) * dt_step, 0.0)
            state['water'] = np.maximum(state['water'] + 0.5 * (d_water + d_water_c) * dt_step, 0.0)
            
            # Acumular masa evacuada de biovolátiles, GNC y vapor
            d_vol_step = np.sum(0.5 * (rates['d_volatiles_dt'] + rates_c['d_volatiles_dt']) * self.dV) * dt_step
            d_gnc_step = np.sum(0.5 * (rates['d_gnc_dt'] + rates_c['d_gnc_dt']) * self.dV) * dt_step
            d_vap_step = np.sum(0.5 * (rates['d_vap_dt'] + rates_c['d_vap_dt']) * self.dV) * dt_step
            
            cum_volatiles_produced += d_vol_step
            cum_gnc_produced += d_gnc_step
            cum_water_evaporated += d_vap_step
            
            # Snapshots
            if next_snap_idx < len(snapshot_times) and t >= snapshot_times[next_snap_idx]:
                snapshots.append({
                    'time': t,
                    'T_map': np.copy(state['T']),
                    'char_map': np.copy(state['char']),
                    'label': f"t = {t/60.0:.1f} min"
                })
                next_snap_idx += 1
                
            # Registro periódico
            if (t - last_save_t) >= save_interval_s or t >= t_final:
                last_save_t = t
                time_history.append(t)
                T_center_history.append(state['T'][0, mid_z_idx])
                T_mid_history.append(state['T'][mid_r_idx, mid_z_idx])
                T_wall_history.append(state['T'][-1, mid_z_idx])
                T_furnace_history.append(T_furnace)
                T_setpoint_history.append(T_set)
                
                M_char = self.integrate_mass(state['char'])
                M_solid = M_char + self.integrate_mass(state['biomass']) + self.integrate_mass(state['water'])
                residual_mass_fraction.append((M_solid / M_init_total) * 100.0)
                
                yield_biochar_history.append((M_char / M_init_total) * 100.0)
                yield_volatiles_history.append((cum_volatiles_produced / M_init_total) * 100.0)
                yield_gnc_history.append((cum_gnc_produced / M_init_total) * 100.0)
                yield_water_history.append((cum_water_evaporated / M_init_total) * 100.0)
                
                # Promedios de tasas volumétricas (kg/(m^3·h))
                rate_biomass_decay_history.append(np.mean(rates['r_biomass_decay']) * 3600.0)
                rate_volatiles_history.append(np.mean(rates['r_vol_prim']) * 3600.0)
                rate_char_history.append(np.mean(rates['r_char_prim']) * 3600.0)
                rate_gnc_history.append(np.mean(rates['r_gnc_prim']) * 3600.0)
                
                if verbose and (len(time_history) % 15 == 0 or t >= t_final):
                    print(
                        f"t: {t/60.0:5.1f} min | "
                        f"T_horno: {T_furnace-273.15:5.1f}°C | "
                        f"T_pared: {state['T'][-1, mid_z_idx]-273.15:5.1f}°C | "
                        f"T_centro: {state['T'][0, mid_z_idx]-273.15:5.1f}°C | "
                        f"Biochar: {yield_biochar_history[-1]:5.1f}% | "
                        f"Volátiles: {yield_volatiles_history[-1]:5.1f}% | "
                        f"GNC: {yield_gnc_history[-1]:5.1f}%"
                    )

            t += dt_step

        if len(snapshots) < len(snapshot_times):
            snapshots.append({
                'time': t_final,
                'T_map': np.copy(state['T']),
                'char_map': np.copy(state['char']),
                'label': f"Final (t = {t_final/60.0:.1f} min)"
            })

        M_unconverted = self.integrate_mass(state['biomass'])
        unconverted_pct = (M_unconverted / M_init_total) * 100.0

        final_yields = {
            'biochar_pct': yield_biochar_history[-1],
            'biovolatiles_pct': yield_volatiles_history[-1],
            'gnc_pct': yield_gnc_history[-1],
            'water_vapor_pct': yield_water_history[-1],
            'unconverted_pct': unconverted_pct
        }

        balance_check = sum(final_yields.values())
        if verbose:
            print(f"-----------------------------------------------------------")
            print(f"RENDIMIENTOS FINALES (MODELO LUMPED):")
            print(f"  • Biochar (sólido carbonoso):              {final_yields['biochar_pct']:.2f} %")
            print(f"  • Biovolátiles / Bio-oil (condensables):   {final_yields['biovolatiles_pct']:.2f} %")
            print(f"  • Gases No Condensables (GNC):             {final_yields['gnc_pct']:.2f} %")
            print(f"  • Vapor de agua (humedad):                 {final_yields['water_vapor_pct']:.2f} %")
            if unconverted_pct > 0.1:
                print(f"  • Biomasa no convertida:                   {unconverted_pct:.2f} %")
            print(f"  TOTAL BALANCE DE MASA:                     {balance_check:.2f} %")
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
            'yield_biovolatiles': np.array(yield_volatiles_history),
            'yield_gnc': np.array(yield_gnc_history),
            'yield_water': np.array(yield_water_history),
            'rate_biomass_decay': np.array(rate_biomass_decay_history),
            'rate_volatiles': np.array(rate_volatiles_history),
            'rate_char': np.array(rate_char_history),
            'rate_gnc': np.array(rate_gnc_history),
            'final_yields': final_yields,
            'snapshots': snapshots,
            'geo': self.geo,
            'ramp': self.ramp
        }
        
        return results
