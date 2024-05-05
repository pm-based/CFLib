function [optionPrice, priceCI] = AsianOptPriceMC(spotPrice, strike, rate, TTM, ...
    putFlag, priceModel, modelParams, nSims, nMonitoring, flagAV)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Add the directory containing the price models to the MATLAB path
addpath(genpath(fullfile('..', 'PriceProcesses')));

if nargin < 10
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
        [X_t, X_t_AV] = MertonProcess(modelParams.sigmaD, modelParams.lambda, ...
            modelParams.muJ, modelParams.sigmaJ, TTM, nMonitoring, nSims, flagAV);
    case 'Kou'
        % Simulate asset prices at maturity using the Kou jump-diffusion model
        [X_t, X_t_AV] = KouProcess(modelParams.sigmaD, modelParams.lambda, ...
            modelParams.lambdaP, modelParams.lambdaN, ...
            modelParams.p, TTM, nMonitoring, nSims, flagAV);
end

% Extract the simulated asset prices at maturity
discountedPayoff = exp(-rate*TTM) * max(0, mean(putFlag * (spotPrice * exp(X_t) - strike), 2));


if flagAV
    % Compute the antithetic payoffs
    discountedPayoff_AV = exp(-rate*TTM) * max(0, mean(putFlag * (spotPrice * exp(X_t_AV) - strike), 2));

    % Calculate the option price and confidence interval using both original and antithetic payoffs
    [optionPrice, ~, priceCI] = normfit((discountedPayoff + discountedPayoff_AV) / 2);
else
    % Calculate the option price and confidence interval using only the original payoffs
    [optionPrice, ~, priceCI] = normfit(discountedPayoff);
end
end

