function [optionPrice, priceCI] = UEOptPriceMC(spotPrice, strike, rate, TTM, ...
    putFlag, priceModel, modelParams, nSims)
% UEOptPriceMC Calculates the price of a European option using Monte Carlo 
% simulation.
%
% This function simulates the final asset price using either the Merton or 
% Kou jump-diffusion model and calculates the European option price and a 
% confidence% interval for the price estimate.
%
% Inputs:
%   spotPrice   - Current price of the underlying asset.
%   strike      - Strike price of the option.
%   rate        - Annual risk-free interest rate, expressed as a decimal.
%   TTM         - Time to maturity of the option, in years.
%   putFlag     - Boolean flag, true for a put option, false for a call option.
%   priceModel  - String specifying the asset price model ('Merton' or 'Kou').
%   modelParams - Struct containing parameters for the chosen price model.
%   nSims       - Number of Monte Carlo simulations to run.
%
% Outputs:
%   optionPrice - Estimated price of the European option.
%   priceCI     - 95% confidence interval for the option price estimate.
%
% The function first adds the 'PriceProcesses' directory to the MATLAB path
% to access necessary pricing models. It then adjusts the putFlag for use 
% in payout calculation, simulates the asset prices at maturity according 
% to the specified model, and finally calculates the option price and 
% confidence interval based on the discounted payoffs.

% Add the directory containing the price models to the MATLAB path
addpath(genpath(fullfile('..', 'PriceProcesses')));

% Convert the putFlag into a numerical value that will be used in the payoff calculation
if (putFlag == true)
    putFlag = -1;  % For a put option, payoff involves max(K - S, 0)
else
    putFlag = 1;   % For a call option, payoff involves max(S - K, 0)
end

% Select and run the appropriate asset price model based on the input
switch priceModel
    case 'Merton'
        % Simulate asset prices at maturity using the Merton jump-diffusion model
        X_T = SimMertonAssetPrice(modelParams.sigmaD, modelParams.lambda, ...
            modelParams.muJ, modelParams.sigmaJ, TTM, 1, nSims);
    case 'Kou'
        % Simulate asset prices at maturity using the Kou jump-diffusion model
        X_T = SimKouAssetPrice(modelParams.sigmaD, modelParams.lambda, ...
                            modelParams.lambdaP, modelParams.lambdaN, ...
                            modelParams.p, TTM, 1, nSims);
end

% Extract the last simulation result as the final asset prices
X_T = X_T(:,end);

% Calculate the European option price and confidence interval using
% the normal approximation of the mean of discounted payoffs
[optionPrice, ~, priceCI] = normfit( exp(-rate*TTM) * max(0, putFlag*(spotPrice*exp(X_T) - strike)) );
end
