function dxdt = funciones_lumped(time, x)
%% FUNCIONES_LUMPED.M
% Sistema de EDOs discretizado en 2D (r, z) para el MODELO CINÉTICO AGRUPADO (Lumped).
% Especies de estado (4 * N variables):
% 1 a N:        Temperatura T (K)
% N+1 a 2N:     Biomasa virgen de café rho_B (kg/m^3)
% 2N+1 a 3N:    Biochar rho_C (kg/m^3)
% 3N+1 a 4N:    Humedad líquida rho_w (kg/m^3)

param_lumped; % Carga parámetros agrupados

dxdt = zeros(4 * N, 1);

T       = x(1 : N);
rho_bio = x(N + 1 : 2 * N);
rho_char= x(2 * N + 1 : 3 * N);
rho_w   = x(3 * N + 1 : 4 * N);

% 1. Horno con inercia de piedra refractaria
if time <= 0
    T_set = T_ini_K;
elseif time <= t_ramp_end
    T_set = T_ini_K + beta_K_s * time;
else
    T_set = T_fin_K;
end
T_furnace = T_ini_K + (T_set - T_ini_K) * (1.0 - exp(-time / tau_stone));

% 2. Propiedades variables y conversión
rho_solid = rho_bio + rho_char + rho_w;
conv = max(0.0, min(1.0, 1.0 - (rho_bio + rho_w) ./ (rho_biomass_0 + rho_water_0)));

k_eff = (1.0 - conv) .* (k_bm0 + k_bm_slope .* (T - 298.15)) + ...
        conv .* (k_char0 + k_char_slope .* (T - 298.15));
k_eff = max(k_eff, 0.05);

cp_eff = (1.0 - conv) .* (1250.0 + 1.8 .* (T - 273.15)) + ...
         conv .* (720.0 + 2.2 .* (T - 273.15));
cp_eff = max(cp_eff, 500.0);
rho_eff = max(rho_solid, 50.0);

% 3. Cinética Lumped Arrhenius
T_clip = max(T, 273.15);
RT = Ru .* T_clip;

k1 = A1 .* exp(-Ea1 ./ RT); % Biomasa -> Biochar
k2 = A2 .* exp(-Ea2 ./ RT); % Biomasa -> Biovolátiles / Bio-oil
k3 = A3 .* exp(-Ea3 ./ RT); % Biomasa -> Gases No Condensables (GNC)
k4 = A4 .* exp(-Ea4 ./ RT); % Volátiles -> Biochar secundario
k5 = A5 .* exp(-Ea5 ./ RT); % Volátiles -> GNC secundario

k_evap = zeros(N, 1);
idx_evap = T_clip > 350.0;
k_evap(idx_evap) = A_evap .* exp(-Ea_evap ./ RT(idx_evap));

% Tasas de reacción (kg/(m^3·s))
rho_B_pos = max(rho_bio, 0.0);
r1 = k1 .* rho_B_pos;
r2 = k2 .* rho_B_pos;
r3 = k3 .* rho_B_pos;
r_evap = k_evap .* max(rho_w, 0.0);

% Craqueo modulado por purga de N2
frac_crk = max(0.0, min(0.30, k5 .* tau_res_vapor));
r5_gnc_sec = frac_crk .* r2;
r4_char_sec = 0.08 .* r5_gnc_sec;

% Derivadas de componentes
d_biomass = -(r1 + r2 + r3);
d_char    = r1 + r4_char_sec;
d_water   = -r_evap;

% Calor de reacción neto (W/m^3)
q_rxn = r_evap .* dH_evap + r1 .* dH1 + r2 .* dH2 + r3 .* dH3 + r4_char_sec .* dH4 + r5_gnc_sec .* dH5;

% 4. Transferencia de calor 2D (r, z)
dT_dt = zeros(N, 1);
dr2 = dr^2;
dz2 = dz^2;

for j = 1 : NL
    for i = 1 : Nr
        k = i + (j - 1) * Nr;
        
        % Axial z
        if j == 1
            d2T_dz2 = (T(k + Nr) - T(k)) / dz2;
        elseif j == NL
            d2T_dz2 = (T(k - Nr) - T(k)) / dz2;
        else
            d2T_dz2 = (T(k + Nr) - 2.0 * T(k) + T(k - Nr)) / dz2;
        end
        term_axial = k_eff(k) * d2T_dz2;
        
        % Radial r
        if i == 1
            % Centro r = 0 (Regla de L'Hôpital): 4*k*(T_1 - T_0)/dr^2
            term_radial = 4.0 * k_eff(k) * (T(k + 1) - T(k)) / dr2;
        elseif i == Nr
            % Pared r = R (Robin con horno)
            q_in = h_comb * (T_furnace - T(k));
            d2T_dr2 = 2.0 * (T(k - 1) - T(k) + (dr * q_in / k_eff(k))) / dr2;
            term_radial = k_eff(k) * (d2T_dr2 + (1.0 / R) * (q_in / k_eff(k)));
        else
            % Nodos internos
            r_curr = r_vec(i);
            d2T_dr2 = (T(k + 1) - 2.0 * T(k) + T(k - 1)) / dr2;
            dT_dr   = (T(k + 1) - T(k - 1)) / (2.0 * dr);
            dk_dr   = (k_eff(k + 1) - k_eff(k - 1)) / (2.0 * dr);
            term_radial = k_eff(k) * (d2T_dr2 + (1.0 / r_curr) * dT_dr) + dk_dr * dT_dr;
        end
        
        dT_dt(k) = (term_radial + term_axial + q_rxn(k)) / (rho_eff(k) * cp_eff(k));
    end
end

% 5. Asignar derivadas
dxdt(1 : N)             = dT_dt;
dxdt(N + 1 : 2 * N)     = d_biomass;
dxdt(2 * N + 1 : 3 * N) = d_char;
dxdt(3 * N + 1 : 4 * N) = d_water;

end
