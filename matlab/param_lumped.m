%% PARAM_LUMPED.M
% Parámetros del modelo cinético agrupado (Lumped Model) para MATLAB.
% Describe la pirólisis de residuos de café mediante 4 especies agrupadas:
% 1. Biomasa virgen de café
% 2. Biochar
% 3. Biovolátiles / Bio-oil (condensables)
% 4. Gases No Condensables (GNC: CO, CO2, CH4, H2)
% 5. Humedad (evaporación a vapor H2O)

%% 1. GEOMETRÍA DEL LECHO Y MALLA
R = 0.025;             % Radio interno (m) [25 mm]
L = 0.20;              % Longitud de la canastilla (m) [20 cm]
Nr = 11;               % Nodos radiales (r)
NL = 15;               % Nodos axiales (z)
N = Nr * NL;           % Total nodos

dr = R / (Nr - 1);
dz = L / (NL - 1);

r_vec = linspace(0, R, Nr);
z_vec = linspace(0, L, NL);

%% 2. CARACTERIZACIÓN DEL CAFÉ (GIEM + FTIR)
rho_bulk = 380.0;      % Densidad aparente empacada (kg/m^3)
w_H2O = 0.08;          % Humedad (8% acondicionado)
w_ASH = 0.0231 * (1 - w_H2O); % Cenizas iniciales (2.31% base seca)

rho_water_0   = rho_bulk * w_H2O;
rho_char_0    = rho_bulk * w_ASH;
rho_biomass_0 = rho_bulk - rho_water_0 - rho_char_0;

%% 3. PARÁMETROS CINÉTICOS AGRUPADOS (LUMPED ARRHENIUS)
Ru = 8.31446;          % J/(mol·K)

% Biomasa -> Biochar primario
A1 = 5.6e7;            Ea1 = 106.5e3;       dH1 = 120.0e3;  % J/kg

% Biomasa -> Biovolátiles / Bio-oil (condensables)
A2 = 1.1e10;           Ea2 = 148.0e3;       dH2 = -210.0e3; % J/kg

% Biomasa -> Gases No Condensables primarios (GNC)
A3 = 4.4e7;            Ea3 = 111.0e3;       dH3 = -80.0e3;  % J/kg

% Biovolátiles -> Biochar secundario (repolimerización)
A4 = 1.0e5;            Ea4 = 95.0e3;        dH4 = 60.0e3;   % J/kg

% Biovolátiles -> Gases No Condensables (craqueo térmico)
A5 = 4.28e6;           Ea5 = 108.0e3;       dH5 = -42.0e3;  % J/kg

% Evaporación de agua
A_evap = 5.0e5;        Ea_evap = 41.8e3;    dH_evap = -2.26e6; % J/kg

tau_res_vapor = 3.0;   % Tiempo de residencia con purga de N2 (s)

%% 4. TERMOFÍSICA Y HORNO
k_bm0 = 0.12;          k_bm_slope = 1.5e-4;
k_char0 = 0.08;        k_char_slope = 2.0e-4;
h_comb = 85.0;         % W/(m^2·K)
tau_stone = 180.0;     % Inercia térmica de la piedra refractaria (s)

%% 5. PROGRAMA DE RAMPA TÉRMICA
beta_C_min = 10.0;     % Tasa de calentamiento (°C/min)
T_ini_C    = 25.0;     % Temperatura inicial (°C)
T_fin_C    = 500.0;    % Temperatura máxima (°C)
t_soak_min = 30.0;     % Tiempo de sostenimiento (min)

T_ini_K = T_ini_C + 273.15;
T_fin_K = T_fin_C + 273.15;
beta_K_s = beta_C_min / 60.0;

t_ramp_end = (T_fin_K - T_ini_K) / beta_K_s;
t_final = t_ramp_end + (t_soak_min * 60.0);
