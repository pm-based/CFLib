function [optionPrice, priceCI] = UEOptPriceMC(spotPrice, strike, rate, TTM, ...
    putFlag, priceModel, modelParams, nSims, VarReduction, barrierType, barrier, nMonitoring)
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

if nargin < 10      % If barrier type is not defined
    barrier = nan;
    barrierType = 'none';
    nMonitoring = 1;
    barrierCheck = @(S_t,barrier) 1;

elseif nargin < 12 % If barrier or nMonitoring value is not defined
    barrier = nan;
    barrierType = 'none';
    nMonitoring = 1;
    barrierCheck = @(S_t,barrier) 1;
    disp('ERROR: no barrier or nMonitoring value specified, no barrier used.')

elseif barrierType == 'DO'
    barrierCheck = @(S_t,barrier) (min(S_t,[],2) > barrier);
    
elseif barrierType == 'DI'
    barrierCheck = @(S_t,barrier) (min(S_t,[],2) < barrier);

elseif barrierType == 'UO'
    barrierCheck = @(S_t,barrier) (max(S_t,[],2) > barrier);

elseif barrierType == 'UI'
    barrierCheck = @(S_t,barrier) (max(S_t,[],2) < barrier);

else
    barrier = nan;
    barrierType = 'none';
    barrierCheck = @(S_t,barrier) 1;
    disp('ERROR: barrier type not allowed, no barrier used.')
end

if nargin < 9       % If the neither VarReduction nor the others arg are defined
    % Variance Reduction args
    VarReduction = 'none';
    flagAV = false;

elseif ~strcmp(VarReduction, 'AV')
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
discountedPayoff = exp(-rate*TTM) * max(0, putFlag * (S_t(:,end) - strike)) .* barrierCheck(S_t,barrier);


switch VarReduction
    case 'AV'
        % Compute the antithetic payoffs
        discountedPayoff_AV = exp(-rate*TTM) * max(0, putFlag * (S_t_AV(:,end) - strike)) .* barrierCheck(S_t_AV,barrier);

        % Calculate the option price and confidence interval using both original and antithetic payoffs
        [optionPrice, ~, priceCI] = normfit((discountedPayoff + discountedPayoff_AV) / 2);

    case 'CV'
        % 1. sample alpha
        nSims_CV = nSims/100;
        S_t_CV = SimAssetPrice(spotPrice, rate, TTM, nMonitoring, nSims_CV, priceModel, modelParams, false);
        f = S_t_CV(:,end);
        g = exp(-rate*TTM) * max(putFlag * (f-strike), 0) .* barrierCheck(S_t_CV,barrier);
        VC = cov(f,g);
        alpha = -VC(1,2)/VC(1,1);

        % 2. compute the price
        f = S_t(:,end);
        Ef = spotPrice * exp(rate*TTM);
        g = exp(-rate*TTM) * max(putFlag * (f-strike), 0) .* barrierCheck(S_t,barrier);
        [optionPrice,~,priceCI] = normfit( g + alpha*(f-Ef) );

    otherwise
        % Calculate the option price and confidence interval using only the original payoffs
        [optionPrice, ~, priceCI] = normfit(discountedPayoff);

end
end
