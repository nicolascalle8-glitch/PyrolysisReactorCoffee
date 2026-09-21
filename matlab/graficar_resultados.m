%% GRAFICAR_RESULTADOS.M
% Generador de gráficos para MATLAB:
% 1. Mapa de calor 2D axisimétrico (contourf) del lecho cilíndrico
% 2. Perfiles dinámicos de temperatura en el tiempo
% 3. Curva de masa residual (% peso vs tiempo)
% 4. Gráfico de barras de rendimientos finales

%% 1. MAPA DE CALOR 2D AL FINAL DE LA PIRÓLISIS
figure('Name', 'Mapa de Calor 2D del Reactor', 'Color', 'w', 'Position', [100, 100, 600, 700]);
T_final_matrix = reshape(C(end, 1 : N), [Nr, NL]) - 273.15; % en °C

[R_grid_mm, Z_grid_cm] = meshgrid(r_vec * 1000.0, z_vec * 100.0);
contourf(R_grid_mm, Z_grid_cm, T_final_matrix', 30, 'LineColor', 'none');
colormap('inferno');
c = colorbar;
ylabel(c, 'Temperatura (°C)', 'FontSize', 11, 'FontWeight', 'bold');
hold on;
[C_lines, h_lines] = contour(R_grid_mm, Z_grid_cm, T_final_matrix', 8, 'LineColor', 'w', 'LineWidth', 0.8);
clabel(C_lines, h_lines, 'FontSize', 9, 'Color', 'w');

xlabel('Radio r (mm)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Altura axial z (cm)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Mapa Térmico 2D en Lecho de Café (t = %.1f min)', t_final/60), 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%% 2. CURVAS DE TEMPERATURA EN EL TIEMPO
figure('Name', 'Evolución Térmica', 'Color', 'w', 'Position', [150, 150, 750, 500]);
t_min = t / 60.0;
mid_z = round(NL / 2);

% Nodos a media altura
k_center = 1 + (mid_z - 1) * Nr;
k_mid    = round(Nr / 2) + (mid_z - 1) * Nr;
k_wall   = Nr + (mid_z - 1) * Nr;

T_center = C(:, k_center) - 273.15;
T_mid    = C(:, k_mid) - 273.15;
T_wall   = C(:, k_wall) - 273.15;

% Setpoint teórico
T_set_curve = zeros(size(t));
for s = 1 : length(t)
    if t(s) <= t_ramp_end
        T_set_curve(s) = T_ini_C + beta_C_min * (t(s) / 60.0);
    else
        T_set_curve(s) = T_fin_C;
    end
end

plot(t_min, T_set_curve, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Setpoint Programado (PID)');
hold on;
plot(t_min, T_wall, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Pared exterior (r = R)');
plot(t_min, T_mid, 'g-.', 'LineWidth', 1.8, 'DisplayName', 'Radio medio (r = R/2)');
plot(t_min, T_center, 'b-', 'LineWidth', 2.2, 'DisplayName', 'Centro (r = 0, L''Hôpital)');

xlabel('Tiempo (minutos)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Temperatura (°C)', 'FontSize', 11, 'FontWeight', 'bold');
title('Perfiles de Temperatura en el Lecho Vertical de Café', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'SouthEast', 'FontSize', 10);
grid on;

%% 3. CURVA TERMOGRAVIMÉTRICA DE PÉRDIDA DE MASA
figure('Name', 'Pérdida de Masa Residual', 'Color', 'w', 'Position', [200, 200, 700, 450]);
M_fraction = (M_solid_t ./ M_init_total) * 100.0;
plot(t_min, M_fraction, 'k-', 'LineWidth', 2.5);
xlabel('Tiempo (minutos)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Masa residual (% wt)', 'FontSize', 11, 'FontWeight', 'bold');
title('Curva de Pérdida de Masa (Termogravimetría en Horno)', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0, 105]);
grid on;

%% 4. DIAGRAMA DE BARRAS DE RENDIMIENTOS FINALES
figure('Name', 'Rendimientos Finales', 'Color', 'w', 'Position', [250, 250, 650, 450]);
categories = {'Biochar', 'Bio-oil / Volátiles', 'Biogás', 'Vapor H2O'};
values = [Yield_Biochar, Yield_Biooil, Yield_Biogas, Yield_Water];

b = bar(values, 0.55, 'FaceColor', [0.2, 0.4, 0.6], 'EdgeColor', 'k', 'LineWidth', 1.2);
set(gca, 'XTickLabel', categories, 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Rendimiento (% en peso)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Distribución de Productos de Pirólisis (Rampa %.0f °C/min)', beta_C_min), 'FontSize', 12, 'FontWeight', 'bold');
ylim([0, max(values) * 1.25]);
grid on;

% Añadir etiquetas con el valor encima de cada barra
xtips = b.XEndPoints;
ytips = b.YEndPoints;
labels = string(round(values, 1)) + "%";
text(xtips, ytips, labels, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 11, 'FontWeight', 'bold');
