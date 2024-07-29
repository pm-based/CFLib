function [S_t, t_i, S_t_AV] = SimAssetPrice(spotPrice, rate, T, nTimeSteps, nProcesses, priceModel, modelParams, flagAV)
%SIMASSETPRICE Simulate asset prices using specified jump-diffusion models
% This function simulates the prices of an asset over time, given initial conditions 
% and parameters for selected jump-diffusion models. It can utilize different models
% as specified by the user, returning the simulated asset prices and the corresponding time steps.
%
% INPUTS:
%   spotPrice   - Initial price of the asset at t = 0.
%   rate        - Risk-free interest rate, used in the exponential adjustment.
%   T           - Total time to maturity of the asset price simulation.
%   nTimeSteps  - Number of time steps to divide the total time T into.
%   nProcesses  - Number of simulated paths or processes.
%   priceModel  - Name of the jump-diffusion model to use ('Merton' or 'Kou').
%   modelParams - Structure containing model-specific parameters such as:
%                 sigmaD (volatility of the diffusion component),
%                 lambda (intensity rate of the jump component),
%                 muJ, sigmaJ (parameters for the jump size distribution in Merton),
%                 lambdaP, lambdaN, p (parameters for the jump size distribution in Kou).
%   flagAV      - Boolean flag indicating whether to compute antithetic variates (true/false).
%
% OUTPUTS:
%   S_t - Simulated asset prices at each time step across all processes.
%   t_i - Vector of time points at which asset prices are simulated.
%
% EXAMPLES:
%   modelParams = struct('sigmaD', 0.1, 'lambda', 0.2, 'muJ', -0.05, 'sigmaJ', 0.1);
%   [S_t, t_i] = SimAssetPrice(100, 0.05, 1, 100, 1000, 'Merton', modelParams);

if nargin < 8
    flagAV = false;
end

switch priceModel
    case 'GBM'
        [X_t, t_i, X_t_AV] = GBMProcess(modelParams.sigma, T, nTimeSteps, nProcesses, flagAV);
    case 'Merton'
        % Simulate asset prices using the Merton jump-diffusion model.
        [X_t, t_i, X_t_AV] = MertonProcess(modelParams.sigmaD, modelParams.lambda, ...
            modelParams.muJ, modelParams.sigmaJ, T, nTimeSteps, nProcesses, flagAV);
    case 'Kou'
        % Simulate asset prices using the Kou jump-diffusion model.
        [X_t, t_i, X_t_AV] = KouProcess(modelParams.sigmaD, modelParams.lambda, ...
                            modelParams.lambdaP, modelParams.lambdaN, ...
                            modelParams.p, T, nTimeSteps, nProcesses, flagAV);
end

% Compute the asset prices at each time step using the formula:
% S_t = S_0 * exp(-r * t_i + X_t), where S_0 is the initial spot price, 
% r is the risk-free rate, t_i are the time points, and X_t is the process from the selected model.
S_t    = spotPrice * exp(rate * t_i + X_t);
if flagAV == true
    S_t_AV = spotPrice * exp(rate * t_i + X_t_AV);
else
    S_t_AV = [];
end
