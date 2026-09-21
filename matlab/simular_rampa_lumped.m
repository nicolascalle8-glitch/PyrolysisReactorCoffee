%% SIMULAR_RAMPA_LUMPED.M
% Script principal para ejecutar la simulación con el MODELO CINÉTICO AGRUPADO (Lumped).
% Resuelve el sistema 2D mediante ode15s de MATLAB.

clear; clc; close all;
tic;

fprintf('===========================================================\n');
fprintf('  SIMULADOR DE PIRÓLISIS: MODELO CINÉTICO AGRUPADO (LUMPED)\n');
fprintf('===========================================================\n');

% 1. Cargar parámetros agrupados
param_lumped;

fprintf('Geometría: R = %.1f mm, L = %.1f cm (Nodos: Nr=%d, NL=%d)\n', R*1000, L*100, Nr, NL);
fprintf('Rampa programada: %.1f °C/min hasta %.1f °C (Soak: %.1f min)\n', beta_C_min, T_fin_C, t_soak_min);
fprintf('Tiempo total: %.1f min (%.0f s)\n', t_final/60, t_final);

% 2. Vector de condiciones iniciales (4 * N variables)
C0 = zeros(4 * N, 1);
C0(1 : N)             = T_ini_K;
C0(N + 1 : 2 * N)     = rho_biomass_0;
C0(2 * N + 1 : 3 * N) = rho_char_0;
C0(3 * N + 1 : 4 * N) = rho_water_0;

% 3. Integración con ode15s
opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'MaxStep', 5.0);
fprintf('Resolviendo con ode15s...\n');
[t, C] = ode15s(@funciones_lumped, [0, t_final], C0, opts);
fprintf('Integración completada en %.2f segundos.\n', toc);

% 4. Postprocesamiento de rendimientos (Yields)
% Volúmenes de control diferenciales para coordenadas cilíndricas dV = 2*pi*r*dr*dz
dV = zeros(Nr, NL);
for i = 1 : Nr
    if i == 1
        area_ring = pi * (dr / 2.0)^2;
    elseif i == Nr
        area_ring = pi * (R^2 - (R - dr/2.0)^2);
    else
        area_ring = 2.0 * pi * r_vec(i) * dr;
    end
    dV(i, :) = area_ring * dz;
end
dV_vec = dV(:);

M_init_total = rho_bulk * pi * R^2 * L;

num_steps = length(t);
M_char_t    = zeros(num_steps, 1);
M_solid_t   = zeros(num_steps, 1);
M_biomass_t = zeros(num_steps, 1);
M_water_t   = zeros(num_steps, 1);

for s = 1 : num_steps
    bio_s  = C(s, N + 1 : 2 * N)';
    char_s = C(s, 2 * N + 1 : 3 * N)';
    w_s    = C(s, 3 * N + 1 : 4 * N)';
    
    M_biomass_t(s) = sum(bio_s .* dV_vec);
    M_char_t(s)    = sum(char_s .* dV_vec);
    M_water_t(s)   = sum(w_s .* dV_vec);
    M_solid_t(s)   = sum((bio_s + char_s + w_s) .* dV_vec);
end

% Rendimientos finales
Yield_Biochar     = (M_char_t(end) / M_init_total) * 100.0;
Yield_Water       = (w_H2O) * 100.0;
Yield_Unconverted = (M_biomass_t(end) / M_init_total) * 100.0;

% Fracción de biomasa consumida que se convierte en Volátiles y Gases No Condensables (GNC)
M_bio_consumed = (M_init_total * (1 - w_H2O - w_ASH)) - M_biomass_t(end);
% Relación estequiométrica entre volátiles y GNC según k2 y k3 a temperatura media de pirólisis (~450°C)
ratio_vol = 0.465;
ratio_gnc = 0.535;

M_volatile_total = M_bio_consumed * (1 - 0.32) * ratio_vol; % Biovolátiles netos
M_gnc_total      = M_bio_consumed * (1 - 0.32) * ratio_gnc; % Gases No Condensables

Yield_Biovolatiles = (M_volatile_total / M_init_total) * 100.0;
Yield_GNC          = max(0.0, 100.0 - Yield_Biochar - Yield_Biovolatiles - Yield_Water - Yield_Unconverted);

fprintf('\n-----------------------------------------------------------\n');
fprintf('RENDIMIENTOS FINALES (MODELO LUMPED):\n');
fprintf('  • Biochar (sólido carbonoso):              %5.2f %%\n', Yield_Biochar);
fprintf('  • Biovolátiles / Bio-oil (condensables):   %5.2f %%\n', Yield_Biovolatiles);
fprintf('  • Gases No Condensables (GNC):             %5.2f %%\n', Yield_GNC);
fprintf('  • Vapor de agua (humedad inicial):         %5.2f %%\n', Yield_Water);
if Yield_Unconverted > 0.1
    fprintf('  • Biomasa no convertida:                   %5.2f %%\n', Yield_Unconverted);
end
fprintf('  BALANCE TOTAL DE MASA:                     %5.2f %%\n', Yield_Biochar + Yield_Biovolatiles + Yield_GNC + Yield_Water + Yield_Unconverted);
fprintf('===========================================================\n');

% 5. Graficar
graficar_lumped;
