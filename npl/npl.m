function [Ak, Pnpl] = npl(G, N0, gamma, S, N, M, Delta)
%NPL  Implements a Nonlinear Power Loading (NPL) algorithm for resource allocation.
%   AUTHOR: [JAKUB KASJANOWICZ,JULIUSZ BOJARCZUK,GRZEGORZ STEPNIAK / Warsaw University of Technology]
%
%   [AK, PNPL] = NPL(G, N0, GAMMA, S, N, M, DELTA) performs a power-loading
%   procedure where power is iteratively allocated across N subcarriers.
%   The algorithm limits total power to N*M (Plim) by incrementally
%   distributing small amounts (Delta) of power to whichever subcarrier
%   yields the maximum capacity increment.
%
%   Inputs:
%       G       - Channel linear response.
%       N0      - Noise spectral density.
%       gamma   - LED nonlinearity coefficient.
%       S       - Nonlinear distortion power spectral distribution.
%       N       - Number of subcarriers.
%       M       - Average power multiplier.
%       Delta   - Power increment step size.
%
%   Outputs:
%       AK      - Final power-loading amplitudes across subcarriers.
%       PNPL    - Final power-loading values across subcarriers.
%
%   The algorithm ends when the total allocated power reaches (or exceeds)
%   the power limit Plim. Capacity is estimated by summing
%   log2(1 + SNRproj), and power is allocated to the subcarrier with the
%   highest incremental capacity gain.
%
%   See also: WATERFILL

%% ------------------ Input Validation ------------------

% Make sure G is a numeric vector.
if ~isvector(G) || ~isnumeric(G)
    error('Input G must be a numeric vector.');
end
if ~isvector(S) || ~isnumeric(S)
    error('S must be a numeric vector.');
end
% Convert G and S to a column vector for consistency.
G = G(:)';
S = S(:)';

% Validate N0, gamma, M, and Delta as numeric scalars.
if ~isscalar(N0) || ~isnumeric(N0)
    error('N0 must be a numeric scalar.');
end
if ~isscalar(gamma) || ~isnumeric(gamma)
    error('gamma must be a numeric scalar.');
end
if ~isscalar(M) || ~isnumeric(M)
    error('M must be a numeric scalar.');
end
if ~isscalar(Delta) || ~isnumeric(Delta)
    error('Delta must be a numeric scalar.');
end

% Validate N as a positive integer matching the length of G.
if ~isscalar(N) || ~isnumeric(N) || (N ~= floor(N)) || (N <= 0)
    error('N must be a positive integer.');
end
if length(G) ~= N || length(S) ~= N
    error('Length of G and S must match N.');
end

%% ------------------ Algorithm Setup ------------------

% Calculate the power limit based on number of subcarriers and a per-subcarrier limit
Plim = N * M;

% Current power distribution across subcarriers (start with zero)
PCurr = zeros(1, N);

% Initialize the final NPL power distribution
Pnpl = zeros(1, N);

% Vector for storing capacity (C) increments for each subcarrier
C = zeros(1, N);

% Square of channel gains for convenience in repeated calculations
G2 = G.^2;

% Iteration counter (optional, can be removed or used for debugging)
it = 1;

%% ------------------ Main Loop ------------------
% Repeat until the total allocated power reaches or exceeds the limit
while sum(PCurr) < Plim
    % Keep a snapshot of the current power distribution to store as final in each step
    Pnpl = PCurr;

    % Evaluate the incremental capacity for distributing Delta power to each subcarrier
    for n = 1 : N
        % Copy the current power distribution
        Pproj = PCurr;

        % Increment power on the nth subcarrier
        Pproj(n) = Pproj(n) + Delta;

        % Compute the projected SNR for all subcarriers given this new power distribution
        SNRproj = Pproj .* G2 ./ (N0 + gamma * S * (sum(Pproj .* G2)).^2);

        % Compute the total capacity (sum of log2(1 + SNR)) with this projection
        C(n) = sum(log2(1 + SNRproj));
    end

    % Find which subcarrier (nn) offers the highest capacity increment
    [~, nn] = max(C);

    % Allocate an additional Delta power to that subcarrier
    PCurr(nn) = PCurr(nn) + Delta;

    % Increment iteration counter
    it = it + 1;
end

%% ------------------ Final Outputs ------------------
% Final allocated amplitudes (square root of the final power)
Ak = sqrt(Pnpl);

end

