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

%
dt = T/nTimeSteps;
X_t = zeros(nProcesses,nTimeSteps+1);

Psi = @(u) -0.5*(sigmaD*u)^2 + 1i*u*lambda*(p/(lambdaP-1i*u) - (1-p)/(lambdaN+1i*u));
muD = -Psi(-1i);

% Diffusion part
Z = randn(nProcesses,nTimeSteps);

% Jumps part
Ndt = icdf('Poisson',rand(nProcesses,nTimeSteps),lambda*dt);

for j = 1:nTimeSteps
    X_t(:,j+1) = X_t(:,j) + muD*dt + sigmaD * sqrt(dt) * Z(:,j);
    for i = 1:nProcesses
        if (Ndt(i,j) > 0)
            Y = 0;
            nPositiveJumps = sum( rand(Ndt(i,j),1) < p );
            nNegativeJumps = Ndt(i,j) - nPositiveJumps;
            if (nPositiveJumps > 0)
                Y = sum(icdf('Exponential',rand(nPositiveJumps,1),1/lambdaP));
            elseif (nNegativeJumps > 0) 
                Y = Y - sum(icdf('Exponential',rand(nNegativeJumps,1),1/lambdaN));
            end
            X_t(i,j+1) = X_t(i,j+1) + Y;
        end
    end
end

