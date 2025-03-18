%% Example Usage of the NPL (Nonlinear Power Loading) Function
addpath("npl\");addpath("results\");

% Set OFDM bandwidth to 100 MHz
BW = 100e6;

% Number of subcarriers
N = 128;

% 3 dB cut-off frequency of the LED
f3dB = 45e6;

% Noise power spectral density
N0 = 1e-4;

% Frequency vector
f = (0:(N-1))/N*BW;

% -------------------------------------------------------------------------
% Channel response modeling
% -------------------------------------------------------------------------
% Here, we define a simple channel response, G, as the magnitude of a 
% low-pass filter: G(f) = 1 / (1 + j*(f / f3dB)).
G = abs(1./(1 + 1i * f / f3dB));

% -------------------------------------------------------------------------
% Nonlinearity modeling
% -------------------------------------------------------------------------
% LED nonlinearity coefficient - describes how strong nonlinearity is
gamma = 1e-11;

% fA is an arbitrary frequency parameter used in this example
fA = 1e6;
% S is a nonlinear distortion power spectral distribution
% (see: 10.1109/JLT.2021.3129586)
% We normalize it by its maximum value.
S = (f3dB^4 * (4*f.^2 + fA^2)) ./ (fA^2 * (4*f.^2 + f3dB^2) .* (f.^2 + f3dB^2));
S = S / max(S);

% -------------------------------------------------------------------------
% NPL parameters
% -------------------------------------------------------------------------
% M is an average power scaling factor.
M = 40;

% Power increment step.
Delta = 0.5;

% -------------------------------------------------------------------------
% Call the NPL algorithm
% -------------------------------------------------------------------------
[Ak, Pnpl] = npl(G, N0, gamma, S, N, M, Delta);

% -------------------------------------------------------------------------
% Plot Results
% -------------------------------------------------------------------------
% Plot the final power allocation across subcarriers
f1 = figure('color','w');
plot(1:N,Pnpl,'o');
xlim([1 N])
title(['N=' num2str(N) ', M=' num2str(M), ', \gamma=', num2str(gamma)])
xlabel('Subcarrier number [#]')
ylabel('Allocated power [a.u.]')
exportgraphics(f1,'results/PowerAlloc.png')