function [optionPrice, priceCI] = UEOptPriceMC(spotPrice, strike, rate, TTM, ...
    putFlag, priceModel, modelParams, nSims, flagAV)
% UEOptPriceMC Calculates the price of a European option using Monte Carlo simulation.
%
% This function simulates the final asset price using either the Merton or Kou
% jump-diffusion model and calculates the European option price and a confidence
% interval for the price estimate. It optionally supports the use of antithetic
% variates to reduce variance.
%
% INPUTS:
%   spotPrice   - Current price of the underlying asset.
%   strike      - Strike price of the option.
%   rate        - Annual risk-free interest rate, expressed as a decimal.
%   TTM         - Time to maturity of the option, in years.
%   putFlag     - Boolean flag, true for a put option, false for a call option.
%   priceModel  - String specifying the asset price model ('Merton' or 'Kou').
%   modelParams - Struct containing parameters for the chosen price model.
%   nSims       - Number of Monte Carlo simulations to run.
%   flagAV      - Boolean flag indicating whether to compute antithetic variates (true/false).
%
% OUTPUTS:
%   optionPrice - Estimated price of the European option.
%   priceCI     - 95% confidence interval for the option price estimate.
%
% DESCRIPTION:
%   The function first adds the 'PriceProcesses' directory to the MATLAB path
%   to access necessary pricing models. It then adjusts the putFlag for use 
%   in payout calculation, simulates the asset prices at maturity according 
%   to the specified model, and finally calculates the option price and 
%   confidence interval based on the discounted payoffs. The optional antithetic
%   variates feature is used to reduce variance in Monte Carlo simulations.
%
% EXAMPLES:
%   % Example with Merton model
%   modelParams = struct('sigmaD', 0.2, 'muJ', -0.1, 'sigmaJ', 0.1, 'lambda', 0.2);
%   [optionPrice, priceCI] = UEOptPriceMC(100, 105, 0.05, 1, false, 'Merton', modelParams, 10000, true);
%   % This example estimates the price of a call option using 10,000 simulations
%   % with antithetic variates enabled.
%
%   % Example with Kou model
%   modelParams = struct('sigmaD', 0.2, 'lambda', 0.1, 'lambdaP', 10, 'lambdaN', 20, 'p', 0.6);
%   [optionPrice, priceCI] = UEOptPriceMC(100, 105, 0.05, 1, true, 'Kou', modelParams, 10000, false);
%   % This example estimates the price of a put option using 10,000 simulations
%   % without antithetic variates.

% Add the directory containing the price models to the MATLAB path
addpath(genpath(fullfile('..', 'PriceProcesses')));

if nargin < 9
    flagAV = false;
end

% Convert the putFlag into a numerical value that will be used in the payoff calculation
if putFlag
    putFlag = -1;  % For a put option, payoff involves max(K - S, 0)
else
    putFlag = 1;   % For a call option, payoff involves max(S - K, 0)
end

% Select and run the appropriate asset price model based on the input
switch priceModel
    case 'Merton'
        % Simulate asset prices at maturity using the Merton jump-diffusion model
        [X_T, X_T_AV] = MertonProcess(modelParams.sigmaD, modelParams.lambda, ...
            modelParams.muJ, modelParams.sigmaJ, TTM, 1, nSims, flagAV);
    case 'Kou'
        % Simulate asset prices at maturity using the Kou jump-diffusion model
        [X_T, X_T_AV] = KouProcess(modelParams.sigmaD, modelParams.lambda, ...
                            modelParams.lambdaP, modelParams.lambdaN, ...
                            modelParams.p, TTM, 1, nSims, flagAV);
end

% Extract the simulated asset prices at maturity
X_T = X_T(:, end);
discountedPayoff = exp(-rate*TTM) * max(0, putFlag * (spotPrice * exp(X_T) - strike));

if flagAV
    % Compute the antithetic payoffs
    X_T_AV = X_T_AV(:, end);
    discountedPayoff_AV = exp(-rate*TTM) * max(0, putFlag * (spotPrice * exp(X_T_AV) - strike));
    
    % Calculate the option price and confidence interval using both original and antithetic payoffs
    [optionPrice, ~, priceCI] = normfit((discountedPayoff + discountedPayoff_AV) / 2);
else
    % Calculate the option price and confidence interval using only the original payoffs
    [optionPrice, ~, priceCI] = normfit(discountedPayoff);
end

end
