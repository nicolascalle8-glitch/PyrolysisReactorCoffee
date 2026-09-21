function dxdt = funciones_ct(time, x)
%% FUNCIONES_CT.M
% Sistema de EDOs discretizado en 2D (r, z) para el reactor vertical de pirólisis.
% Resuelve simultáneamente:
% 1. Ecuación de energía con L'Hôpital en r = 0 y convección en r = R
% 2. Balances de masa para agua, celulosa, hemicelulosa, lignina, lípidos y biochar.
% Compatible con ode15s de MATLAB.

param; % Carga parámetros geométricos, cinéticos y termofísicos

dxdt = zeros(8 * N, 1);

% Extracción de variables de estado
T        = x(1 : N);
rho_w    = x(N + 1 : 2 * N);
rho_cell = x(2 * N + 1 : 3 * N);
rho_hce  = x(3 * N + 1 : 4 * N);
rho_lig  = x(4 * N + 1 : 5 * N);
rho_lip  = x(5 * N + 1 : 6 * N);
rho_char = x(6 * N + 1 : 7 * N);
rho_tar  = x(7 * N + 1 : 8 * N);

% -------------------------------------------------------------------------
% 1. TEMPERATURA DEL HORNO CON INERCIA DE PIEDRA REFRACTARIA
% -------------------------------------------------------------------------
if time <= 0
    T_set = T_ini_K;
elseif time <= t_ramp_end
    T_set = T_ini_K + beta_K_s * time;
else
    T_set = T_fin_K;
end
% Dinámica de primer orden de la piedra refractaria:
T_furnace = T_ini_K + (T_set - T_ini_K) * (1.0 - exp(-time / tau_stone));

% -------------------------------------------------------------------------
% 2. PROPIEDADES TERMOFÍSICAS VARIABLES
% -------------------------------------------------------------------------
rho_solid = rho_w + rho_cell + rho_hce + rho_lig + rho_lip + rho_char;
conv = max(0.0, min(1.0, 1.0 - rho_solid ./ rho_bulk));

k_eff = (1.0 - conv) .* (k_bm0 + k_bm_slope .* (T - 298.15)) + ...
        conv .* (k_char0 + k_char_slope .* (T - 298.15));
k_eff = max(k_eff, 0.05);

cp_eff = (1.0 - conv) .* (1250.0 + 1.8 .* (T - 273.15)) + ...
         conv .* (720.0 + 2.2 .* (T - 273.15));
cp_eff = max(cp_eff, 500.0);

rho_eff = max(rho_solid, 50.0);

% -------------------------------------------------------------------------
% 3. CINÉTICA QUÍMICA DE REACCIÓN
% -------------------------------------------------------------------------
T_clip = max(T, 250.0);
RT = Ru .* T_clip;

% Velocidades de reacción de Arrhenius (1/s)
k_evap = zeros(N, 1);
idx_evap = T_clip > 350.0;
k_evap(idx_evap) = A_evap .* exp(-Ea_evap ./ RT(idx_evap));

k_hce  = A_hce  .* exp(-Ea_hce  ./ RT);
k_cell = A_cell .* exp(-Ea_cell ./ RT);
k_lig  = A_lig  .* exp(-Ea_lig  ./ RT);
k_lip  = A_lip  .* exp(-Ea_lip  ./ RT);
k_crk  = A_crk_g.* exp(-Ea_crk_g./ RT);

% Tasas de consumo (kg/(m^3·s))
r_evap = k_evap .* max(rho_w, 0.0);
r_hce  = k_hce  .* max(rho_hce, 0.0);
r_cell = k_cell .* max(rho_cell, 0.0);
r_lig  = k_lig  .* max(rho_lig, 0.0);
r_lip  = k_lip  .* max(rho_lip, 0.0);

% Craqueo limitado por arrastre de gas inerte N2 (tiempo de residencia tau ~ 3 s)
frac_crk = max(0.0, min(0.35, k_crk .* 3.0));
r_crk_g = frac_crk .* (y_tar_hce.*r_hce + y_tar_cell.*r_cell + y_tar_lig.*r_lig + y_tar_lip.*r_lip);
r_crk_c = 0.08 .* r_crk_g;

% Producción neta de biochar y volátiles
d_char = y_char_hce.*r_hce + y_char_cell.*r_cell + y_char_lig.*r_lig + y_char_lip.*r_lip + r_crk_c;
d_tar  = max(0.0, y_tar_hce.*r_hce + y_tar_cell.*r_cell + y_tar_lig.*r_lig + y_tar_lip.*r_lip - r_crk_g - r_crk_c);

% Término fuente de calor por reacción volumétrica (W/m^3)
q_rxn = r_evap .* dH_evap + r_hce .* dH_hce + r_cell .* dH_cell + r_lig .* dH_lig + r_lip .* dH_lip;

% -------------------------------------------------------------------------
% 4. ECUACIÓN DE TRANSFERENCIA DE CALOR 2D DISCRETIZADA
% -------------------------------------------------------------------------
dT_dt = zeros(N, 1);
dr2 = dr^2;
dz2 = dz^2;

for j = 1 : NL
    for i = 1 : Nr
        k = i + (j - 1) * Nr; % Índice del nodo actual
        
        % Conducción axial (z)
        if j == 1
            d2T_dz2 = (T(k + Nr) - T(k)) / dz2; % Frontera superior
        elseif j == NL
            d2T_dz2 = (T(k - Nr) - T(k)) / dz2; % Frontera inferior
        else
            d2T_dz2 = (T(k + Nr) - 2.0 * T(k) + T(k - Nr)) / dz2;
        end
        term_axial = k_eff(k) * d2T_dz2;
        
        % Conducción radial (r)
        if i == 1
            % =============================================================
            % CENTRO (r = 0): REGLA DE L'HÔPITAL (PÁGINA 3 DOCUMENTO)
            % lim_{r->0} (1/r)*dT/dr = d2T/dr2  => 4*k*(T_1 - T_0)/dr^2
            % =============================================================
            term_radial = 4.0 * k_eff(k) * (T(k + 1) - T(k)) / dr2;
            
        elseif i == Nr
            % =============================================================
            % PARED EXTERIOR (r = R): CONVECCIÓN Y RADIACIÓN DEL HORNO
            % -k * dT/dr = h_comb * (T_wall - T_furnace)
            % =============================================================
            q_in = h_comb * (T_furnace - T(k));
            d2T_dr2 = 2.0 * (T(k - 1) - T(k) + (dr * q_in / k_eff(k))) / dr2;
            term_radial = k_eff(k) * (d2T_dr2 + (1.0 / R) * (q_in / k_eff(k)));
            
        else
            % =============================================================
            % NODOS INTERNOS (1 < i < Nr)
            % =============================================================
            r_curr = r_vec(i);
            d2T_dr2 = (T(k + 1) - 2.0 * T(k) + T(k - 1)) / dr2;
            dT_dr   = (T(k + 1) - T(k - 1)) / (2.0 * dr);
            dk_dr   = (k_eff(k + 1) - k_eff(k - 1)) / (2.0 * dr);
            
            term_radial = k_eff(k) * (d2T_dr2 + (1.0 / r_curr) * dT_dr) + dk_dr * dT_dr;
        end
        
        % Balance de energía en el nodo
        dT_dt(k) = (term_radial + term_axial + q_rxn(k)) / (rho_eff(k) * cp_eff(k));
    end
end

% -------------------------------------------------------------------------
% 5. ASIGNACIÓN DEL VECTOR DE DERIVADAS
% -------------------------------------------------------------------------
dxdt(1 : N)             = dT_dt;
dxdt(N + 1 : 2 * N)     = -r_evap;
dxdt(2 * N + 1 : 3 * N) = -r_cell;
dxdt(3 * N + 1 : 4 * N) = -r_hce;
dxdt(4 * N + 1 : 5 * N) = -r_lig;
dxdt(5 * N + 1 : 6 * N) = -r_lip;
dxdt(6 * N + 1 : 7 * N) = d_char;
dxdt(7 * N + 1 : 8 * N) = d_tar;

end
