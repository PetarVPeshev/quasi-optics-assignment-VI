close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../quasi-optics-library');
c = physconst('LightSpeed');

%% PARAMETERS
wave.f = linspace(6, 14, 1001) * 1e9;
medium.er = 1;
dipole.l = 14e-3;
dipole.w = 1e-3;
fss.dx = 15e-3;
fss.dy = 15e-3;
theta = 0;
phi = 0;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k = 2 * pi ./ wave.wavelength;

%% SPHERICAL COORDINATE GRID
sph_grid = NaN(1, 1, 2);
sph_grid(:, :, 1) = theta;
sph_grid(:, :, 2) = phi;

Gamma = NaN(1, length(wave.f));
T = NaN(1, length(wave.f));
for idx = 1 : 1 : length(wave.f)
    %% WAVE VECTOR COMPONENTS
    [k_comp, ~] = wave_vector(medium.er, wave.k(idx), sph_grid);

    %% BF CURRENT
    ibf = bf_current(medium.er, dipole.l, dipole.w, fss.dx, fss.dy, ...
        wave.k(idx), sph_grid);

    %% REFLECTION AND TRANSMISSION COEFFICIENTS
    [Gamma(idx), T(idx)] = fss_reflection(medium.er, dipole.l, ...
        dipole.w, fss.dx, fss.dy, ibf, wave.k(idx), k_comp, sph_grid);
end

%% PLOT REFLECTION AND TRANSMISSION COEFFICIENTS
figure('Position', [250 250 800 400]);
plot(wave.f * 1e-9, abs(Gamma), 'LineWidth', 2.0, 'DisplayName', '\Gamma');
hold on;
plot(wave.f * 1e-9, abs(T), '--', 'LineWidth', 2.0, 'DisplayName', 'T');
hold on;
plot(wave.f * 1e-9, (abs(Gamma) .^ 2) + (abs(T) .^ 2), ...
    'LineWidth', 2.0, 'DisplayName', '|\Gamma|^{2} + |T|^{2}');
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('f / GHz');
ylabel('|\Gamma|, |T|');
