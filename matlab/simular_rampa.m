%% SIMULAR_RAMPA.M
% Script principal para ejecutar la simulación de transferencia de calor
% y pirólisis de residuos de café (SCG) en MATLAB.
% Resuelve el sistema 2D mediante el integrador implícito ode15s.

clear; clc; close all;
tic;

fprintf('===========================================================\n');
fprintf('  SIMULADOR DE PIRÓLISIS EN HORNO VERTICAL (MATLAB)\n');
fprintf('===========================================================\n');

% 1. Cargar parámetros del sistema
param;

fprintf('Geometría: R = %.1f mm, L = %.1f cm (Nodos: Nr=%d, NL=%d)\n', R*1000, L*100, Nr, NL);
fprintf('Rampa programada: %.1f °C/min hasta %.1f °C (Soak: %.1f min)\n', beta_C_min, T_fin_C, t_soak_min);
fprintf('Tiempo total de simulación: %.1f min (%.0f s)\n', t_final/60, t_final);
fprintf('Iniciando integración con ode15s...\n');

% 2. Construir vector de condiciones iniciales (tamaño 8*N)
C0 = zeros(8 * N, 1);
C0(1 : N)             = T_ini_K;
C0(N + 1 : 2 * N)     = rho_w_0;
C0(2 * N + 1 : 3 * N) = rho_cell_0;
C0(3 * N + 1 : 4 * N) = rho_hce_0;
C0(4 * N + 1 : 5 * N) = rho_lig_0;
C0(5 * N + 1 : 6 * N) = rho_lip_0;
C0(6 * N + 1 : 7 * N) = rho_char_0;
C0(7 * N + 1 : 8 * N) = 0.0; % Volátiles iniciales

% 3. Configuración y resolución de EDOs rígidas (ode15s)
opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'MaxStep', 5.0);
[t, C] = ode15s(@funciones_ct, [0, t_final], C0, opts);

fprintf('Integración finalizada en %.2f segundos.\n', toc);

% 4. Postprocesamiento y cálculo de rendimientos (Yields)
% Volúmenes de control diferenciales para integración cilíndrica dV = 2*pi*r*dr*dz
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
M_char_t = zeros(num_steps, 1);
M_tar_t  = zeros(num_steps, 1);
M_solid_t= zeros(num_steps, 1);
M_w_t    = zeros(num_steps, 1);

for s = 1 : num_steps
    w_s    = C(s, N + 1 : 2 * N)';
    cell_s = C(s, 2 * N + 1 : 3 * N)';
    hce_s  = C(s, 3 * N + 1 : 4 * N)';
    lig_s  = C(s, 4 * N + 1 : 5 * N)';
    lip_s  = C(s, 5 * N + 1 : 6 * N)';
    char_s = C(s, 6 * N + 1 : 7 * N)';
    tar_s  = C(s, 7 * N + 1 : 8 * N)';
    
    M_w_t(s)    = sum(w_s .* dV_vec);
    M_char_t(s) = sum(char_s .* dV_vec);
    M_tar_t(s)  = sum(tar_s .* dV_vec);
    M_solid_t(s)= sum((w_s + cell_s + hce_s + lig_s + lip_s + char_s) .* dV_vec);
end

% Rendimientos finales (% en peso)
Yield_Biochar = (M_char_t(end) / M_init_total) * 100.0;
Yield_Biooil  = (M_tar_t(end) / M_init_total) * 100.0;
Yield_Water   = (w_H2O) * 100.0; % Evaporación completa
Yield_Biogas  = max(0.0, 100.0 - Yield_Biochar - Yield_Biooil - Yield_Water);

fprintf('\n-----------------------------------------------------------\n');
fprintf('RENDIMIENTOS FINALES ESTIMADOS:\n');
fprintf('  • Biochar (sólido carbonoso):         %5.2f %%\n', Yield_Biochar);
fprintf('  • Biovolátiles / Bio-oil (alquitranes):%5.2f %%\n', Yield_Biooil);
fprintf('  • Biogás (gases no condensables):     %5.2f %%\n', Yield_Biogas);
fprintf('  • Vapor de agua (humedad inicial):    %5.2f %%\n', Yield_Water);
fprintf('  BALANCE TOTAL DE MASA:                %5.2f %%\n', Yield_Biochar + Yield_Biooil + Yield_Biogas + Yield_Water);
fprintf('===========================================================\n');

% 5. Generar gráficas
graficar_resultados;
