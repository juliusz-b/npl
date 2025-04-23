function [Ak, Pnpl] = npl_opt(G, N0, gamma, S, N, M)
%NPL_OPT  Implements Nonlinear Power Loading using optimization techniques.
%   AUTHOR: [JULIUSZ BOJARCZUK / Warsaw University of Technology]
%
%   [AK, PNPL] = NPL_OPT(G, N0, GAMMA, S, N, M) performs power-loading
%   using MATLAB's optimization functions to maximize capacity subject to 
%   exact power constraints. 
%
%   Inputs:
%       G       - Channel linear response.
%       N0      - Noise spectral density.
%       gamma   - LED nonlinearity coefficient.
%       S       - Nonlinear distortion power spectral distribution.
%       N       - Number of subcarriers.
%       M       - Average power multiplier.
%
%   Outputs:
%       AK      - Final power-loading amplitudes across subcarriers.
%       PNPL    - Final power-loading values across subcarriers.
%
%   The algorithm directly optimizes the power allocation across all
%   subcarriers simultaneously.
%
%   See also: NPL, WATERFILL, FMINCON

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

% Validate N0, gamma, and M as numeric scalars.
if ~isscalar(N0) || ~isnumeric(N0)
    error('N0 must be a numeric scalar.');
end
if ~isscalar(gamma) || ~isnumeric(gamma)
    error('gamma must be a numeric scalar.');
end
if ~isscalar(M) || ~isnumeric(M)
    error('M must be a numeric scalar.');
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

% Square of channel gains for convenience in calculations
G2 = G.^2;

% Initial guess - equal power with slight variation to avoid local minima
P0 = ones(1, N) * (M * 0.8) + rand(1, N) * (M * 0.2);
P0 = P0 * (Plim / sum(P0));

% Define optimization constraints
lb = zeros(1, N);           % Lower bound: power must be non-negative
ub = ones(1, N) * Plim;     % Upper bound: no single carrier can exceed total power
Aeq = ones(1, N);           % Equality constraint matrix: sum of all powers
beq = Plim;                 % Equality constraint value: must equal Plim
A = [];                     % No inequality constraints 
b = [];                     % No inequality constraints

% Optimization options for fmincon with gradient specification
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'off', ...
    'MaxFunctionEvaluations', 1e6, ...
    'MaxIterations', 2000, ...
    'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-14, ...
    'ConstraintTolerance', 1e-12, ...
    'SpecifyObjectiveGradient', true);  % Enable analytical gradient use

%% ------------------ Optimization Execution ------------------

% Optimize power allocation to maximize capacity with gradient information
[Pnpl, ~, exitflag] = fmincon(@(P) neg_capacity_with_grad(P, G2, N0, gamma, S), ...
    P0, A, b, Aeq, beq, lb, ub, [], options);

% Check if optimization was successful
if exitflag <= 0
    warning('Optimization might not have converged. Exit flag: %d', exitflag);
end

%% ------------------ Final Outputs ------------------

% Final allocated amplitudes (square root of the final power)
Ak = sqrt(Pnpl);

end

function [neg_cap, grad] = neg_capacity_with_grad(P, G2, N0, gamma, S)
%NEG_CAPACITY_WITH_GRAD  Computes negative capacity and its gradient
%   Calculates the negative channel capacity for a given power allocation
%   and its gradient with respect to each power value.

% Plin
X = sum(P .* G2);

% SNR at all subc.
SNR = P .* G2 ./ (N0 + gamma * S * X^2);

% Capacity
capacity = sum(log2(1 + SNR));
neg_cap = -capacity;

% Grad.
if nargout > 1
    N = length(P);
    grad = zeros(1, N);
    
    distortion_terms = (N0 + gamma * S * X^2);
    
    for k = 1:N
        direct = G2(k) / distortion_terms(k);
        
        indirect = 0;
        for i = 1:N
            dSNR_i = -2 * gamma * S(i) * X * G2(k) * P(i) * G2(i) / (distortion_terms(i)^2);
            indirect = indirect + dSNR_i / ((1 + SNR(i)) * log(2));
        end
        
        grad(k) = -(direct / ((1 + SNR(k)) * log(2)) + indirect);
    end
end
end