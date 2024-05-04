function X_t = KouProcess(sigmaD, lambda, lambdaP, lambdaN, p, T, nTimeSteps, nProcesses)
% KOUPROCESS Simulates paths of a Kou jump-diffusion process.
% This function generates a matrix of simulated values from a Kou jump-diffusion
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
%
% OUTPUTS:
%   X_t          - A matrix of size (nProcesses x (nTimeSteps+1)) where each row represents
%                  a simulated path of the process.

% Calculate the time increment.
dt = T/nTimeSteps;

% Initialize the matrix to store the process values. First column is zero (initial value).
X_t = zeros(nProcesses, nTimeSteps+1);

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
    
    % Loop over each process to apply the jump component.
    for i = 1:nProcesses
        % Check if there are jumps in this time step for the current process.
        if (Ndt(i, j) > 0)
            Y = 0;
            % Determine the number of positive and negative jumps.
            nPositiveJumps = sum(rand(Ndt(i, j), 1) < p);
            nNegativeJumps = Ndt(i, j) - nPositiveJumps;
            
            % Sum the jump sizes using exponentially distributed magnitudes.
            if (nPositiveJumps > 0)
                Y = sum(icdf('Exponential', rand(nPositiveJumps, 1), 1/lambdaP));
            end
            if (nNegativeJumps > 0)
                Y = Y - sum(icdf('Exponential', rand(nNegativeJumps, 1), 1/lambdaN));
            end
            
            % Update the process value by the total jump magnitude.
            X_t(i, j+1) = X_t(i, j+1) + Y;
        end
    end
end
