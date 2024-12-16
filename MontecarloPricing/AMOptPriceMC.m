function [optionPrice, priceCI] = AMOptPriceMC(spotPrice, strike, rate, TTM, ...
    putFlag, priceModel, modelParams, nSims, nMonitoring)
% RIGA 72: 
% TODO: DA ADATTARE ALLA CALL


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

% Convert the putFlag into a numerical value that will be used in the payoff calculation
if putFlag
    putFlag = -1;  % For a put option, payoff involves max(K - S, 0)
else
    putFlag = 1;   % For a call option, payoff involves max(S - K, 0)
end
flagAV = false;

S_t = SimAssetPrice(spotPrice, rate, TTM, nMonitoring, nSims, priceModel, modelParams, flagAV);
S_t=S_t(:,2:end);   % we do not consider t=0 for early exercise

% AM price computation according to Longstuff paper
Exercise_Time = nMonitoring*ones(nSims,1);
payoff = max(0, putFlag * (S_t(:,end) - strike)); %payoff
dt = TTM/nMonitoring;

for j = nMonitoring-1:-1:1
    Inmoney = find(S_t(:,j)<strike);
    S_I = S_t(Inmoney,j);
    %-- Intrinsic Value
    IV = strike-S_I; % DA ADATTARE PER LA CALL
    %-- Continuation Value 
    %- Regression

    % In the regression we need to describe 3 different polynomes:
    % - ones(length(S_I),1) allows to set the blank number
    % - S_I allows to set the first order
    % - S_I.^2 allows to set the second order

    A = [ones(length(S_I),1), S_I, S_I.^2]; 
    idd = Exercise_Time(Inmoney)-j;
    b=payoff(Inmoney).*exp(-rate*dt*(idd));

    % we introduce "Exercise_Time(Inmoney)-j" since we want to look only at
    % those steps when we're in the money from the time j to the time T.
    % Look carefully at the steps during the computation !!
    % It's a count of how many times the k-step is IN THE MONEY !!!

    alpha=A\b;
    %- Continuation Value 
    CV=A*alpha;

    %----------
    %== j is an exercise instant?
    Index=find(IV>CV);
    Early_Exercise=Inmoney(Index);
    % Update
    payoff(Early_Exercise)=IV(Index);
    Exercise_Time(Early_Exercise)=j;
end
[optionPrice, ~ ,priceCI] = normfit(payoff.*exp(-rate*dt*Exercise_Time));

end
