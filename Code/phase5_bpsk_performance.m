%% Phase 5 - BPSK BER performance over simulated fading channels
clear; clc; close all;

%% --- Simulation parameters ---
Nbits = 2e5;                 % number of bits (increase if you want smoother curves)
SNRdB = 0:2:20;              % Eb/N0 in dB
rng(1);                      % repeatable results

% Time axis for channel generation
fs = 1e4;                    % sampling frequency used in your channel simulation
Tsim = 0.2;                  % seconds
t = 0:1/fs:Tsim;

% To map symbols to channel samples: pick random time indices
% (gives i.i.d. samples from the simulated process)
% Alternatively: you can take sequential samples if you want time correlation.

%% --- Channel parameters from your previous steps ---
% Urban (Rayleigh)
fd_urban = 241;              % Hz
sigma_ray_urban = 0.9018;    % Rayleigh sigma

% Highway (Rician)
fd_highway = 578;            % Hz
s_highway = 0.9194;
sigma_ric_highway = 0.3485;

% Sum-of-sinusoids parameters
N1 = 10; N2 = 11;

%% --- Function to generate sum-of-sinusoids Gaussian process (I/Q) ---
gen_gauss_sos = @(fd, sigma, t) local_sos_gauss(fd, sigma, t, N1, N2);

%% --- Generate fading sequences (one realization) ---
% Urban Rayleigh: h = hI + j hQ
[hI_u, hQ_u] = gen_gauss_sos(fd_urban, sigma_ray_urban, t);
h_urban = hI_u + 1j*hQ_u;

% Highway Rician: h = s + (hI + j hQ)
[hI_h, hQ_h] = gen_gauss_sos(fd_highway, sigma_ric_highway, t);
h_highway = s_highway + (hI_h + 1j*hQ_h);

h_urban = h_urban(:);
h_highway = h_highway(:);
%% --- Normalize average channel power to 1 (important for fair SNR) ---
h_urban   = h_urban / sqrt(mean(abs(h_urban).^2));
h_highway = h_highway / sqrt(mean(abs(h_highway).^2));

%% --- Generate bits and BPSK symbols ---
bits = randi([0 1], Nbits, 1);
x = 2*bits - 1;              % BPSK: 0->-1, 1->+1

%% --- Pre-allocate BER arrays ---
BER_urban = zeros(size(SNRdB));
BER_highway = zeros(size(SNRdB));

%% --- Indices into channel samples (random sampling) ---
L = length(t);
idx = randi([1 L], Nbits, 1);
hU = h_urban(idx);
hH = h_highway(idx);

hU = hU(:);
hH = hH(:);
x  = x(:);
bits = bits(:);
%% --- BER loop ---
for ii = 1:length(SNRdB)
    EbN0 = 10^(SNRdB(ii)/10);     % linear
    N0 = 1/EbN0;                  % since Es=Eb=1 for BPSK with unit power
    noise_sigma = sqrt(N0/2);

    % Urban
    nU = noise_sigma*(randn(Nbits,1) + 1j*randn(Nbits,1));
    yU = hU.*x + nU;

    % Coherent equalization (perfect CSI)
    zU = yU ./ hU;
    xhatU = real(zU) > 0;
    BER_urban(ii) = mean(xhatU ~= bits);

    % Highway
    nH = noise_sigma*(randn(Nbits,1) + 1j*randn(Nbits,1));
    yH = hH.*x + nH;

    zH = yH ./ hH;
    xhatH = real(zH) > 0;
    BER_highway(ii) = mean(xhatH ~= bits);

    fprintf('SNR=%2d dB | BER urban=%.3e | BER highway=%.3e\n', ...
        SNRdB(ii), BER_urban(ii), BER_highway(ii));
end

%% --- Plot BER curves ---
figure;
semilogy(SNRdB, BER_urban, 'o-', 'LineWidth', 1.5); hold on; grid on;
semilogy(SNRdB, BER_highway, 's--', 'LineWidth', 1.5);
xlabel('E_b/N_0 (dB)');
ylabel('BER');
title('BPSK Performance over Simulated V2V Fading Channels');
legend('Urban (Rayleigh)', 'Highway (Rician)', 'Location', 'southwest');

%% --- Optional: save figure ---
saveas(gcf, fullfile('..','Figures','ber_curves.png'));

%% ========= Local function =========
function [hI, hQ] = local_sos_gauss(fd, sigma, t, N1, N2)
    % Sum-of-sinusoids Gaussian process generation (Appendix)
    f1 = fd * sin(pi*((1:N1)-0.5)/(2*N1));
    f2 = fd * sin(pi*((1:N2)-0.5)/(2*N2));

    theta1 = 2*pi*rand(1,N1);
    theta2 = 2*pi*rand(1,N2);

    c1 = sigma*sqrt(2/N1);
    c2 = sigma*sqrt(2/N2);

    hI = zeros(size(t));
    hQ = zeros(size(t));

    for n=1:N1
        hI = hI + c1*cos(2*pi*f1(n)*t + theta1(n));
    end
    for n=1:N2
        hQ = hQ + c2*cos(2*pi*f2(n)*t + theta2(n));
    end
end
