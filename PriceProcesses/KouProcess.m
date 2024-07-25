function [X_t, t_i, X_t_AV] = KouProcess(sigmaD, lambda, lambdaP, lambdaN, p, T, nTimeSteps, nProcesses, flagAV)
% KOUPROCESS Simulates paths of a Kou jump-diffusion process.
% This function generates matrices of simulated values from a Kou jump-diffusion
% process, which is useful for modeling financial data with both continuous
% diffusive behaviors and discontinuous jump behaviors.
%
% INPUTS:
%   sigmaD       - Volatility of the diffusion part.
%   lambda       - Average rate of jumps per unit time.
%   lambdaP      - Rate parameter for the exponential distribution of positive jumps.
%   lambdaN      - Rate parameter for the exponential distribution of negative jumps.
%   p            - Probability of a jump being positive.
%   T            - Total time span of the simulation.
%   nTimeSteps   - Number of time steps in the simulation.
%   nProcesses   - Number of independent processes to simulate.
%   flagAV       - Boolean flag indicating whether to compute antithetic variates (true/false).
%
% OUTPUTS:
%   X_t          - A matrix of size (nProcesses x (nTimeSteps+1)) where each row represents
%                  a simulated path of the Kou model.
%   X_t_AV       - A matrix of the same size as X_t, containing the antithetic paths.
%                  If flagAV is false, this will be returned as an empty matrix.
%   t_i          - A vector of length (nTimeSteps+1) representing the discrete time points
%                  at which the process is evaluated.
%
% DESCRIPTION:
%   The Kou jump-diffusion model includes both a continuous Gaussian component (diffusion part)
%   and a jump component characterized by a double exponential distribution. Positive jumps
%   follow one exponential distribution (with rate lambdaP), while negative jumps follow another
%   exponential distribution (with rate lambdaN), making this model asymmetric and suitable for
%   modeling assets with skewed jump behavior. The optional antithetic variates feature is used
%   to reduce variance in Monte Carlo simulations.
%
% EXAMPLES:
%   [X_t, X_t_AV, t_i] = KouProcess(0.2, 0.1, 10, 20, 0.6, 1, 100, 1000, true);
%   % This example generates 1000 simulated paths of a Kou process over 1 year with
%   % 100 time steps, a diffusion volatility of 0.2, a jump intensity of 0.1,
%   % positive and negative jump rates of 10 and 20 respectively, and returns
%   % antithetic paths.

if nargin < 9
    flagAV = false;
end

% Calculate the time increment.
dt = T / nTimeSteps;
t_i = 0:dt:T;

% Initialize the matrices to store the process values. First column is zero (initial value).
X_t = zeros(nProcesses, nTimeSteps+1);
if flagAV
    X_t_AV = zeros(nProcesses, nTimeSteps+1);
else
    X_t_AV = [];
end

% Define the characteristic exponent Psi of the Lévy process associated with the Kou model.
Psi = @(u) -0.5*(sigmaD*u)^2 + 1i*u*lambda*(p/(lambdaP-1i*u) - (1-p)/(lambdaN+1i*u));

% Calculate the drift adjustment to ensure the process is a martingale.
muD = -Psi(-1i);

% Generate random numbers for the Gaussian components (diffusion part).
Z = randn(nProcesses, nTimeSteps);

% Generate Poisson-distributed random counts for the jumps (scaled by dt).
Ndt = icdf('Poisson', rand(nProcesses, nTimeSteps), lambda*dt);

% Loop over each time step to evolve the process.
for j = 1:nTimeSteps
    % Update each process for diffusion and drift.
    X_t(:, j+1) = X_t(:, j) + muD*dt + sigmaD * sqrt(dt) * Z(:, j);

    if flagAV
        X_t_AV(:, j+1) = X_t_AV(:, j) + muD*dt - sigmaD * sqrt(dt) * Z(:, j);
    end

    % Loop over each process to apply the jump component.
    for i = 1:nProcesses
        % Check if there are jumps in this time step for the current process.
        if Ndt(i, j) > 0
            Y = 0;
            Y_AV = 0; % Always generated
            % Determine the number of positive and negative jumps.
            nPositiveJumps = sum(rand(Ndt(i, j), 1) < p);
            nNegativeJumps = Ndt(i, j) - nPositiveJumps;

            % Sum the jump sizes using exponentially distributed magnitudes.
            if nPositiveJumps > 0
                u = rand(nPositiveJumps, 1);
                Y = sum(icdf('Exponential', u, 1/lambdaP));
                if flagAV
                    Y_AV = sum(icdf('Exponential', 1-u, 1/lambdaP));
                end
            end
            if nNegativeJumps > 0
                u = rand(nNegativeJumps, 1);
                Y = Y - sum(icdf('Exponential', u, 1/lambdaN));
                if flagAV
                    Y_AV = Y_AV - sum(icdf('Exponential', 1-u, 1/lambdaN));
                end
            end

            % Update the process value by the total jump magnitude.
            X_t(i, j+1) = X_t(i, j+1) + Y;
            if flagAV
                X_t_AV(i, j+1) = X_t_AV(i, j+1) + Y_AV;
            end
        end
    end
end

end
