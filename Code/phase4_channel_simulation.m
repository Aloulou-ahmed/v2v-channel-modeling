%% Phase 4 - Sum-of-Sinusoids Channel Simulation
clear; clc; close all;

% Parameters
fc = 5.2e9;
fd_urban = 241;      % Hz
fs = 1e4;            % simulation sampling frequency
Tsim = 0.2;          % simulation duration (seconds)
t = 0:1/fs:Tsim;

% Rayleigh sigma (from fitting)
sigma_ray = 0.9018;

% Number of sinusoids
N1 = 10; N2 = 11;

% Frequencies
f1 = fd_urban * sin(pi*( (1:N1)-0.5 )/(2*N1));
f2 = fd_urban * sin(pi*( (1:N2)-0.5 )/(2*N2));

% Random phases
theta1 = 2*pi*rand(1,N1);
theta2 = 2*pi*rand(1,N2);

% Coefficients
c1 = sigma_ray*sqrt(2/N1);
c2 = sigma_ray*sqrt(2/N2);

% Generate I and Q
hI = zeros(size(t));
hQ = zeros(size(t));

for n = 1:N1
    hI = hI + c1*cos(2*pi*f1(n)*t + theta1(n));
end

for n = 1:N2
    hQ = hQ + c2*cos(2*pi*f2(n)*t + theta2(n));
end

% Complex fading channel
h_rayleigh = hI + 1j*hQ;

% Envelope
env_rayleigh = abs(h_rayleigh);

% Plot
figure;
plot(t, env_rayleigh);
xlabel('Time (s)');
ylabel('|h(t)|');
title('Urban Rayleigh Fading Envelope');
grid on;


%% Rician fading - Highway

fd_highway = 578;    % Hz

% Rician parameters (from fitting)
s = 0.9194;
sigma_ric = 0.3485;

% Frequencies
f1 = fd_highway * sin(pi*( (1:N1)-0.5 )/(2*N1));
f2 = fd_highway * sin(pi*( (1:N2)-0.5 )/(2*N2));

theta1 = 2*pi*rand(1,N1);
theta2 = 2*pi*rand(1,N2);

c1 = sigma_ric*sqrt(2/N1);
c2 = sigma_ric*sqrt(2/N2);

hI = zeros(size(t));
hQ = zeros(size(t));

for n = 1:N1
    hI = hI + c1*cos(2*pi*f1(n)*t + theta1(n));
end

for n = 1:N2
    hQ = hQ + c2*cos(2*pi*f2(n)*t + theta2(n));
end

% Add LOS component
h_rician = s + hI + 1j*hQ;

env_rician = abs(h_rician);

figure;
plot(t, env_rician);
xlabel('Time (s)');
ylabel('|h(t)|');
title('Highway Rician Fading Envelope');
grid on;

% Save figures (PNG)
saveas(1, fullfile('..','Figures','rayleigh_fading_envelope.png'));
saveas(2, fullfile('..','Figures','rician_fading_envelope.png'));