%% PARAM.M
% Configuración de parámetros para la simulación del reactor vertical de pirólisis de café
% Geometría, propiedades termofísicas, cinética y rampa térmica

%% 1. PARÁMETROS GEOMÉTRICOS DEL REACTOR
R = 0.025;             % Radio interno del lecho de biomasa (m) [50 mm de diámetro]
L = 0.20;              % Longitud/altura de la canastilla de café (m) [20 cm]
Nr = 11;               % Número de nodos radiales (r)
NL = 15;               % Número de nodos axiales (z)
N = Nr * NL;           % Número total de nodos espaciales

dr = R / (Nr - 1);     % Paso radial (m)
dz = L / (NL - 1);     % Paso axial (m)

% Grilla de posiciones
r_vec = linspace(0, R, Nr);
z_vec = linspace(0, L, NL);

[R_mesh, Z_mesh] = ndgrid(r_vec, z_vec);
r_node = R_mesh(:);    % Vector columna de coordenadas radiales para los N nodos

%% 2. PROPIEDADES DE LA BIOMASA DE RESIDUOS DE CAFÉ (SCG)
% Caracterización obtenida por GIEM (UdeA/EIA) y espectroscopía FTIR:
% - Humedad seca/preparada: 8.0 % wt (fresco = 61.2% según reporte GIEM)
% - Cenizas: 2.31 % wt (base seca)
% - FTIR: Picos de lípidos a 2924 y 2853 cm^-1, éster a 1743 cm^-1, polisacáridos a 1030-1160 cm^-1.
rho_bulk = 380.0;      % Densidad aparente empacada del lecho de café (kg/m^3)
w_H2O = 0.08;          % Fracción de humedad en base húmeda
w_CELL = 0.12 * (1 - w_H2O); % Celulosa
w_HCE  = 0.35 * (1 - w_H2O); % Hemicelulosa
w_LIG  = 0.28 * (1 - w_H2O); % Lignina
w_LIP  = 0.15 * (1 - w_H2O); % Lípidos / Aceites de café
w_ASH  = 0.0231 * (1 - w_H2O); % Cenizas iniciales

% Concentraciones másicas iniciales por metro cúbico de lecho (kg/m^3)
rho_w_0    = rho_bulk * w_H2O;
rho_cell_0 = rho_bulk * w_CELL;
rho_hce_0  = rho_bulk * w_HCE;
rho_lig_0  = rho_bulk * w_LIG;
rho_lip_0  = rho_bulk * w_LIP;
rho_char_0 = rho_bulk * w_ASH;

%% 3. CINÉTICA QUÍMICA DE PIRÓLISIS (Constantes de Arrhenius)
Ru = 8.31446;          % Constante universal de los gases (J/(mol·K))

% Evaporación de agua
A_evap = 5.0e5;        Ea_evap = 41.8e3;     dH_evap = -2.26e6; % J/kg
% Hemicelulosa (HCE)
A_hce  = 1.2e10;       Ea_hce  = 125.0e3;    dH_hce  = -65.0e3;
y_char_hce = 0.28;     y_tar_hce = 0.38;     y_gas_hce = 0.34;
% Celulosa (CELL)
A_cell = 2.5e13;       Ea_cell = 185.0e3;    dH_cell = -150.0e3;
y_char_cell = 0.12;    y_tar_cell = 0.65;    y_gas_cell = 0.23;
% Lignina (LIG)
A_lig  = 8.0e7;        Ea_lig  = 98.0e3;     dH_lig  = 120.0e3;
y_char_lig = 0.52;     y_tar_lig = 0.28;     y_gas_lig = 0.20;
% Lípidos (LIP)
A_lip  = 2.0e9;        Ea_lip  = 120.0e3;    dH_lip  = -180.0e3;
y_char_lip = 0.05;     y_tar_lip = 0.82;     y_gas_lip = 0.13;
% Craqueo de alquitrán
A_crk_g = 4.2e6;       Ea_crk_g = 108.0e3;

%% 4. PROPIEDADES TERMOFÍSICAS Y TRANSPORTE
k_bm0 = 0.12;          % Conductividad biomasa virgen (W/(m·K))
k_bm_slope = 1.5e-4;
k_char0 = 0.08;        % Conductividad biochar (W/(m·K))
k_char_slope = 2.0e-4;

h_comb = 85.0;         % Coeficiente combinado radiación + convección horno-tubo (W/(m^2·K))
rho_gas = 0.8;         % Densidad gas de pirólisis (kg/m^3)
Cp_steam = 2050.0;     % Cp vapor/gases (J/(kg·K))

% Inercia térmica del recubrimiento de piedra del horno (s)
tau_stone = 180.0;

%% 5. PROGRAMA DE RAMPA TÉRMICA POR DEFECTO
% Modifica estas variables según tu ensayo experimental:
beta_C_min = 10.0;     % Tasa de calentamiento (°C/min)
T_ini_C    = 25.0;     % Temperatura inicial (°C)
T_fin_C    = 500.0;    % Temperatura máxima isotérmica (°C)
t_soak_min = 30.0;     % Tiempo de sostenimiento isotérmico (min)

T_ini_K = T_ini_C + 273.15;
T_fin_K = T_fin_C + 273.15;
beta_K_s = beta_C_min / 60.0;

t_ramp_end = (T_fin_K - T_ini_K) / beta_K_s;  % Tiempo de rampa (s)
t_final = t_ramp_end + (t_soak_min * 60.0);   % Tiempo total (s)
