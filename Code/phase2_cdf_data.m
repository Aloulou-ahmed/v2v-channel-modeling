%% Phase 2 - Measured CDF data (Urban + Highway)
clear; clc; close all;

% Levels in dB (same for both scenarios)
levels_dB = [-20 -17.5 -15 -12.5 -10 -7.5 -5 -2.5 0 2.5 5 7.5 10];

% Measured CDF - Urban (Table 2)
cdf_urban = [0.0003 0.0012 0.003 0.007 0.01 0.04 0.06 0.14 0.4 0.8 0.98 0.99 1];

% Measured CDF - Highway (Table 3)
cdf_highway = [0.0004 0.0014 0.0032 0.008 0.02 0.055 0.1 0.225 0.53 0.85 0.95 0.99 1];

%% Plot measured CDFs (sanity check)
figure;
plot(levels_dB, cdf_urban, 'o-', 'LineWidth', 1.5); grid on;
xlabel('Level (dB)');
ylabel('CDF');
title('Measured CDF - Urban Environment');

figure;
plot(levels_dB, cdf_highway, 'o-', 'LineWidth', 1.5); grid on;
xlabel('Level (dB)');
ylabel('CDF');
title('Measured CDF - Highway Environment');

% Save figures (PNG)
saveas(1, fullfile('..','Figures','cdf_urban.png'));
saveas(2, fullfile('..','Figures','cdf_highway.png'));