function [X_t, X_t_AV, t_i] = MertonProcess(sigmaD, muJ, sigmaJ, lambda, T, nTimeSteps, nProcesses, flagAV)
% MERTONPROCESS Simulates paths of a Merton jump-diffusion process.
% This function generates matrices of simulated values from a Merton jump-diffusion
% process, which incorporates both Gaussian diffusion elements and a jump process
% characterized by normally distributed jumps. It optionally supports the use of antithetic variates.
%
% INPUTS:
%   sigmaD       - Volatility of the diffusion part of the stock price.
%   muJ          - Mean of the jump size distribution.
%   sigmaJ       - Standard deviation of the jump size distribution.
%   lambda       - Intensity of the jump process (average number of jumps per unit time).
%   T            - Total time span of the simulation.
%   nTimeSteps   - Number of time steps in the simulation.
%   nProcesses   - Number of independent process paths to simulate.
%   flagAV       - Boolean flag indicating whether to compute antithetic variates (true/false).
%
% OUTPUTS:
%   X_t          - A matrix of size (nProcesses x (nTimeSteps+1)) where each row represents
%                  a simulated path of the Merton model.
%   X_t_AV       - A matrix of the same size as X_t, containing the antithetic paths.
%                  If flagAV is false, this will be returned as an empty matrix.
%   t_i          - A vector of length (nTimeSteps+1) representing the discrete time points
%                  at which the process is evaluated.
%
% DESCRIPTION:
%   The Merton jump-diffusion model includes both a continuous Gaussian diffusion part
%   and a jump component characterized by a Poisson process with normally distributed jump sizes.
%   This allows modeling assets with occasional large jumps in price in addition to the standard
%   continuous diffusive behavior. The optional antithetic variates feature is used to reduce variance
%   in Monte Carlo simulations.

if nargin < 8
    flagAV = false;
end

% Calculate the time increment.
dt = T / nTimeSteps;

% Time vector from 0 to T at intervals of dt.
t_i = dt:dt:T;

% Define the characteristic exponent Psi for the Lévy process in the Merton model.
Psi = @(u) -0.5*(sigmaD*u)^2 ...
            + lambda * (exp(-0.5*(sigmaJ*u)^2 + 1i*muJ*u) - 1);

% Calculate the drift correction to ensure the process is a martingale.
muW = -Psi(-1i);

% Generate random numbers for the Gaussian components (diffusion part).
Z = randn(nProcesses, nTimeSteps);

% Generate Poisson-distributed random counts for the jumps (scaled by dt).
Ndt = icdf('Poisson', rand(nProcesses, nTimeSteps), lambda*dt);

% Compute the cumulative sum of the scaled Gaussian increments to form the Brownian motion.
B_t = sigmaD * cumsum(sqrt(dt) .* Z, 2);

% Initialize a matrix to accumulate jump sizes.
Y = zeros(nProcesses, nTimeSteps);
if flagAV
    Y_AV = zeros(nProcesses, nTimeSteps);
end

% Loop over each process and time step to accumulate jumps.
for i = 1:nProcesses
    for j = 1:nTimeSteps
        if Ndt(i, j) > 0
            % Sum up the jump sizes, where jump magnitudes are normally distributed.
            Y(i, j) = sum(muJ + sigmaJ * randn(1, Ndt(i, j)));
            
            if flagAV
                % Use antithetic variates for the jump sizes.
                Y_AV(i, j) = sum(muJ - sigmaJ * randn(1, Ndt(i, j)));
            end
        end
    end
end

% Combine drift, diffusion, and jump components to compute the process paths.
X_t = muW * t_i .* ones(nProcesses, nTimeSteps) ...
        + B_t ...
        + cumsum(Y, 2);

% Prepend zeros to represent the starting value of the process.
X_t = cat(2, zeros(nProcesses, 1), X_t);
t_i = cat(2, 0, t_i);

if flagAV
    % Compute the antithetic paths with the opposite sign of the diffusion component.
    X_t_AV = muW * t_i .* ones(nProcesses, nTimeSteps) ...
        - B_t ...
        + cumsum(Y_AV, 2);
    
    % Prepend zeros to represent the starting value of the antithetic process.
    X_t_AV = cat(2, zeros(nProcesses, 1), X_t_AV);
else
    % Return an empty matrix if antithetic variates are not requested.
    X_t_AV = [];
end

end
