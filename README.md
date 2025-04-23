# Nonlinear Power Loading (NPL) Algorithm

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15046845.svg)](https://doi.org/10.5281/zenodo.15046845)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=juliusz-b/npl)

## Overview

This repository provides MATLAB implementations of the Nonlinear Power Loading (NPL) algorithm for optical wireless communication systems using LEDs. NPL maximizes transmission capacity while minimizing nonlinear distortion generated in LED transmitters.

Unlike conventional power loading algorithms that treat nonlinear distortion as additive noise, NPL specifically accounts for the nonlinear characteristics of the channel, resulting in up to 10% throughput improvement compared to classical approaches.

## Algorithm Implementations

The repository contains two implementations of the NPL algorithm:

### 1. NPL (Iterative Approach)

```matlab
[Ak, Pnpl] = npl(G, N0, gamma, S, N, M, Delta)
```

This implementation uses an iterative approach where power is incrementally allocated (Delta) to whichever subcarrier yields the maximum capacity increment. The algorithm continues until the total power limit is reached.

### 2. NPL_OPT (Optimization Approach)

```matlab
[Ak, Pnpl] = npl_opt(G, N0, gamma, S, N, M)
```

This implementation uses MATLAB's optimization functions (`fmincon`) with analytical gradients to directly find the optimal power allocation that maximizes capacity. This approach is generally faster and more accurate than the iterative method.

## Parameters

Both functions use the following parameters:

| Parameter | Description |
|-----------|-------------|
| G         | Channel linear response vector |
| N0        | Noise spectral density |
| gamma     | LED nonlinearity coefficient |
| S         | Nonlinear distortion power spectral distribution |
| N         | Number of subcarriers |
| M         | Average power multiplier |
| Delta     | Power increment step size (only for `npl`) |

## Return Values

Both functions return:

| Value | Description |
|-------|-------------|
| Ak    | Final power-loading amplitudes across subcarriers |
| Pnpl  | Final power-loading values across subcarriers |

## Usage Example

```matlab
% Define your system parameters
N = 128;                  % Number of subcarriers
M = 1;                    % Average power multiplier
G = channel_response(N);  % Channel response
N0 = 1e-4;                % Noise spectral density
gamma = 0.1;              % LED nonlinearity coefficient
S = distortion_profile(N);% Distortion power distribution

% For iterative approach
Delta = 0.001;            % Power increment
[Ak_iter, P_iter] = npl(G, N0, gamma, S, N, M, Delta);

% For optimization approach
[Ak_opt, P_opt] = npl_opt(G, N0, gamma, S, N, M);

% Plot and compare results
figure;
subplot(2,1,1);
plot(1:N, P_iter, 'b-', 1:N, P_opt, 'r--');
legend('Iterative', 'Optimization');
title('Power Allocation Comparison');

subplot(2,1,2);
plot(1:N, 10*log10(P_iter./P_opt));
title('Power Difference (dB)');
xlabel('Subcarrier Index');
```

## Citation

TBA

## License

This project is licensed under the LGPL v3 License - see the LICENSE file for details.

## Contributors

- Jakub Kasjanowicz - Warsaw University of Technology
- Juliusz Bojarczuk - Warsaw University of Technology
- Grzegorz Stepniak - Warsaw University of Technology