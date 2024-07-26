function [optionPrice, priceCI] = AsianOptPriceMC(spotPrice, rate, TTM, ...
    putFlag, priceModel, modelParams, nSims, nMonitoring, VarReduction)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Add the directory containing the price models to the MATLAB path
addpath(genpath(fullfile('..', 'PriceProcesses')));

if nargin < 10
    VarReduction = 'none';
    flagAV = false;
elseif VarReduction ~= 'AV'
    flagAV = false;
else
    flagAV = true;
end

% Convert the putFlag into a numerical value that will be used in the payoff calculation
if putFlag
    putFlag = -1;  % For a put option, payoff involves max(K - S, 0)
else
    putFlag = 1;   % For a call option, payoff involves max(S - K, 0)
end

[S_t, ~, S_t_AV] = SimAssetPrice(spotPrice, rate, TTM, nMonitoring, nSims, priceModel, modelParams, flagAV);

% Extract the simulated asset prices at maturity
discountedPayoff = exp(-rate*TTM) * max(0, putFlag *(S_t(:,end) - mean(S_t,2)));


switch VarReduction
    case 'AV'
        % Compute the antithetic payoffs
        discountedPayoff_AV = exp(-rate*TTM) * max(0, putFlag *(S_t_AV(:,end) - mean(S_t_AV,2)));

        % Calculate the option price and confidence interval using both original and antithetic payoffs
        [optionPrice, ~, priceCI] = normfit((discountedPayoff + discountedPayoff_AV) / 2);

    case 'CV'
        % 1. sample alpha
        nSims_CV = nSims/100;
        S_t_CV = SimAssetPrice(spotPrice, rate, TTM, nMonitoring, nSims_CV, priceModel, modelParams, false);
        f = S_t_CV(:,end) - mean(S_t_CV,2);
        g = exp(-rate*TTM) * max(f, 0);
        VC = cov(f,g);
        alpha = -VC(1,2)/VC(1,1);

        % 2. compute the price
        f = S_t(:,end) - mean(S_t,2);
        dt = TTM/nMonitoring;
        Ef = spotPrice * exp(rate*TTM) - mean(spotPrice*exp(rate*(0:nMonitoring)*dt));
        g = exp(-rate*TTM) * max(f, 0);
        [optionPrice,~,priceCI] = normfit( g + alpha*(f-Ef) );

    otherwise
        % Calculate the option price and confidence interval using only the original payoffs
        [optionPrice, ~, priceCI] = normfit(discountedPayoff);

end
end

