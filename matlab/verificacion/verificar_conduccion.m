% verificar_conduccion.m — Prueba 7 del spec del motor nuevo (segundo camino en MATLAB).
%
% Lee caso_conduccion.csv y caso_conduccion_meta.csv (los escribe exportar_caso.py), vuelve a armar
% por su cuenta el sistema de volúmenes finitos (lecho + acero, cilíndricas, caras z adiabáticas,
% resistencia de contacto en serie, Robin en la cara externa), lo integra con ode15s y compara con la
% T que calculó Python (Euler implícito) en los mismos tiempos.
% Criterio: error máximo <= 1 % de (T_horno - T0).
%
% Uso: en MATLAB, cd a matlab/verificacion y escribir  verificar_conduccion
% (también corre en GNU Octave).

clear; clc;
aqui = fileparts(mfilename('fullpath'));
if isempty(aqui), aqui = pwd; end

% ---------- Escalares
fid = fopen(fullfile(aqui, 'caso_conduccion_meta.csv'), 'r');
fgetl(fid);                                   % cabecera
meta = struct();
while true
    linea = fgetl(fid);
    if ~ischar(linea), break; end
    partes = strsplit(linea, ',');
    meta.(strtrim(partes{1})) = str2double(partes{2});
end
fclose(fid);
Nr = meta.Nr; Nl = meta.Nl; Nz = meta.Nz;

% ---------- Nodos (cabecera con los tiempos)
fid = fopen(fullfile(aqui, 'caso_conduccion.csv'), 'r');
cab = strsplit(fgetl(fid), ',');
fclose(fid);
D = dlmread(fullfile(aqui, 'caso_conduccion.csv'), ',', 1, 0);
cols_T = find(strncmp(cab, 'T_', 2));
tiempos = cellfun(@(c) str2double(c(3:end)), cab(cols_T));
i_nodo = D(:, 1) + 1;  k_nodo = D(:, 2) + 1;
idx = sub2ind([Nr, Nz], i_nodo, k_nodo);
r = zeros(Nr, 1); z = zeros(Nz, 1);
r(i_nodo) = D(:, 3); z(k_nodo) = D(:, 4);
k = zeros(Nr, Nz); rc = zeros(Nr, Nz); T0 = zeros(Nr, Nz); V_py = zeros(Nr, Nz);
k(idx) = D(:, 6); rc(idx) = D(:, 7); T0(idx) = D(:, 8); V_py(idx) = D(:, 5);
T_py = cell(numel(tiempos), 1);
for j = 1:numel(tiempos)
    Tj = zeros(Nr, Nz); Tj(idx) = D(:, cols_T(j)); T_py{j} = Tj;
end

% ---------- Geometría de los volúmenes (se rehace, no se copia de Python)
caras = zeros(Nr + 1, 1);
caras(2:Nl) = 0.5 * (r(1:Nl-1) + r(2:Nl));
caras(Nl + 1) = meta.R_lecho;
caras(Nl + 2:Nr) = 0.5 * (r(Nl+1:Nr-1) + r(Nl+2:Nr));
caras(Nr + 1) = meta.R_ext;
dz_k = meta.dz * ones(1, Nz); dz_k([1 Nz]) = meta.dz / 2;
V = pi * (caras(2:end).^2 - caras(1:end-1).^2) * dz_k;
fprintf('Volúmenes: diferencia máxima relativa con Python = %.2e\n', max(abs(V(:) - V_py(:)) ./ V_py(:)));

% ---------- Conductancias (W/K)
armonica = @(a, b) 2 * a .* b ./ (a + b);
G_r = (2 * pi * caras(2:Nr)) * dz_k .* armonica(k(1:end-1, :), k(2:end, :)) ./ diff(r);
R_serie = meta.dr_acero ./ k(Nl + 1, :) + meta.R_contacto_extra;   % nodo del lecho sobre la cara; acero a dr_a
G_r(Nl, :) = 2 * pi * meta.R_lecho * dz_k ./ R_serie;
area_z = pi * (caras(2:end).^2 - caras(1:end-1).^2);
G_z = repmat(area_z, 1, Nz - 1) .* armonica(k(:, 1:end-1), k(:, 2:end)) / meta.dz;
G_ext = 2 * pi * meta.R_ext * dz_k * meta.h_ext;

% ---------- Sistema C dT/dt = M T + c   (orden de MATLAB: columna, n = i + (k-1)·Nr)
N = Nr * Nz; n = reshape(1:N, Nr, Nz);
filas = []; cols = []; vals = [];
a = n(1:end-1, :); b = n(2:end, :);
filas = [filas; a(:); b(:); a(:); b(:)]; cols = [cols; b(:); a(:); a(:); b(:)];
vals = [vals; G_r(:); G_r(:); -G_r(:); -G_r(:)];
a = n(:, 1:end-1); b = n(:, 2:end);
filas = [filas; a(:); b(:); a(:); b(:)]; cols = [cols; b(:); a(:); a(:); b(:)];
vals = [vals; G_z(:); G_z(:); -G_z(:); -G_z(:)];
ext = n(Nr, :)';
filas = [filas; ext]; cols = [cols; ext]; vals = [vals; -G_ext(:)];
M = sparse(filas, cols, vals, N, N);
c = zeros(N, 1); c(ext) = G_ext(:) * meta.T_horno;
Cinv = spdiags(1 ./ (rc(:) .* V(:)), 0, N, N);
J = Cinv * M;  f = @(t, y) J * y + Cinv * c;

opciones = odeset('RelTol', 1e-8, 'AbsTol', 1e-6, 'Jacobian', J);
[~, Y] = ode15s(f, [0, tiempos], T0(:), opciones);

% ---------- Comparación
err = 0;
for j = 1:numel(tiempos)
    err = max(err, max(abs(Y(j + 1, :)' - T_py{j}(:))));
end
dT0 = meta.T_horno - meta.T0;
fprintf('Error máximo MATLAB (ode15s) vs Python (Euler implícito, dt = %g s): %.3f K\n', meta.dt_python, err);
if err <= 0.01 * dT0, veredicto = 'PASA'; else, veredicto = 'NO PASA'; end
fprintf('Límite (1 %% de %.1f K): %.3f K  ->  %s\n', dT0, 0.01 * dT0, veredicto);
