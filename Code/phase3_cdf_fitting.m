%% Phase 3 - CDF fitting
clear; clc; close all;

levels_dB = [-20 -17.5 -15 -12.5 -10 -7.5 -5 -2.5 0 2.5 5 7.5 10];
r = 10.^(levels_dB/20);   % linear envelope amplitude

cdf_urban = [0.0003 0.0012 0.003 0.007 0.01 0.04 0.06 0.14 0.4 0.8 0.98 0.99 1];
cdf_highway = [0.0004 0.0014 0.0032 0.008 0.02 0.055 0.1 0.225 0.53 0.85 0.95 0.99 1];

rayleigh_cdf = @(sigma, r) 1 - exp(-(r.^2)/(2*sigma^2));

% --- Rayleigh fit: Urban ---
err_rayleigh_urban = @(sigma) sum((rayleigh_cdf(sigma, r) - cdf_urban).^2);
sigma_ray_urban = fminsearch(err_rayleigh_urban, 1);

fprintf('Urban Rayleigh sigma = %.4f\n', sigma_ray_urban);


% --- Rayleigh fit: Highway ---
err_rayleigh_highway = @(sigma) sum((rayleigh_cdf(sigma, r) - cdf_highway).^2);
sigma_ray_highway = fminsearch(err_rayleigh_highway, 1);

fprintf('Highway Rayleigh sigma = %.4f\n', sigma_ray_highway);

% --- Plot: Urban (measured vs Rayleigh fit) ---
figure;
plot(levels_dB, cdf_urban, 'o', 'LineWidth', 1.5); hold on; grid on;
plot(levels_dB, rayleigh_cdf(sigma_ray_urban, r), '-', 'LineWidth', 1.8);
xlabel('Level (dB)'); ylabel('CDF');
title('Urban: Measured vs Rayleigh Fit');
legend('Measured', sprintf('Rayleigh (\\sigma=%.3f)', sigma_ray_urban), 'Location', 'best');

% --- Plot: Highway (measured vs Rayleigh fit) ---
figure;
plot(levels_dB, cdf_highway, 'o', 'LineWidth', 1.5); hold on; grid on;
plot(levels_dB, rayleigh_cdf(sigma_ray_highway, r), '-', 'LineWidth', 1.8);
xlabel('Level (dB)'); ylabel('CDF');
title('Highway: Measured vs Rayleigh Fit');
legend('Measured', sprintf('Rayleigh (\\sigma=%.3f)', sigma_ray_highway), 'Location', 'best');

% Save figures (PNG)
saveas(1, fullfile('..','Figures','cdf_urban_fit.png'));
saveas(2, fullfile('..','Figures','cdf_highway_fit.png'));


%% --- Rician CDF using noncentral chi-square (robust) ---
% Rician envelope: parameters s (LOS amplitude), sigma (scatter std)
rician_cdf = @(s, sigma, r) ncx2cdf( (r./sigma).^2, 2, (s./sigma).^2 );

% To enforce positivity, we optimize over log-parameters:
% p(1)=log(s), p(2)=log(sigma)
rician_cdf_from_p = @(p, r) rician_cdf(exp(p(1)), exp(p(2)), r);

%% --- Least-squares fits (Rician) ---
% Urban
err_rician_urban = @(p) sum((rician_cdf_from_p(p, r) - cdf_urban).^2);
p0 = [log(0.5) log(0.9)];   % initial guess [s, sigma]
p_rician_urban = fminsearch(err_rician_urban, p0);

s_urban     = exp(p_rician_urban(1));
sigma_urban = exp(p_rician_urban(2));

% Highway
err_rician_highway = @(p) sum((rician_cdf_from_p(p, r) - cdf_highway).^2);
p0 = [log(1.0) log(0.8)];   % initial guess [s, sigma]
p_rician_highway = fminsearch(err_rician_highway, p0);

s_highway     = exp(p_rician_highway(1));
sigma_highway = exp(p_rician_highway(2));

%% --- Compute K-factor (dimensionless) ---
K_urban   = (s_urban^2)   / (2*sigma_urban^2);
K_highway = (s_highway^2) / (2*sigma_highway^2);

fprintf('\n--- Rician Fit Results ---\n');
fprintf('Urban:   s = %.4f, sigma = %.4f, K = %.4f (%.2f dB)\n', ...
    s_urban, sigma_urban, K_urban, 10*log10(K_urban));
fprintf('Highway: s = %.4f, sigma = %.4f, K = %.4f (%.2f dB)\n', ...
    s_highway, sigma_highway, K_highway, 10*log10(K_highway));

%% --- Plots: Measured vs Rayleigh vs Rician (Urban) ---
figure;
plot(levels_dB, cdf_urban, 'o', 'LineWidth', 1.5); hold on; grid on;
plot(levels_dB, rayleigh_cdf(sigma_ray_urban, r), '-', 'LineWidth', 1.8);
plot(levels_dB, rician_cdf(s_urban, sigma_urban, r), '--', 'LineWidth', 1.8);
xlabel('Level (dB)'); ylabel('CDF');
title('Urban: Measured vs Rayleigh vs Rician');
legend('Measured', ...
       sprintf('Rayleigh (\\sigma=%.3f)', sigma_ray_urban), ...
       sprintf('Rician (s=%.3f, \\sigma=%.3f)', s_urban, sigma_urban), ...
       'Location', 'best');

%% --- Plots: Measured vs Rayleigh vs Rician (Highway) ---
figure;
plot(levels_dB, cdf_highway, 'o', 'LineWidth', 1.5); hold on; grid on;
plot(levels_dB, rayleigh_cdf(sigma_ray_highway, r), '-', 'LineWidth', 1.8);
plot(levels_dB, rician_cdf(s_highway, sigma_highway, r), '--', 'LineWidth', 1.8);
xlabel('Level (dB)'); ylabel('CDF');
title('Highway: Measured vs Rayleigh vs Rician');
legend('Measured', ...
       sprintf('Rayleigh (\\sigma=%.3f)', sigma_ray_highway), ...
       sprintf('Rician (s=%.3f, \\sigma=%.3f)', s_highway, sigma_highway), ...
       'Location', 'best');

% Save figures (PNG)
saveas(3, fullfile('..','Figures','cdf_urban_rayleigh_rician.png'));
saveas(4, fullfile('..','Figures','cdf_highway_rayleigh_rician.png'));