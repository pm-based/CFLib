clear
%clc

%% AM vanilla option pricing with PDE (BS model)
% Market params
spotPrice = 1;
rate = 0.01;
sigma = 0.4;

% option params
strike = 1;
TTM = 1;
putFlag = true;

% PDE params
PDEParams.M = 200; % time 
PDEParams.N = 1000; % logprice


optionPrice = AMOptBSPricePDE(spotPrice, strike, rate, ...
        sigma, TTM, putFlag, PDEParams)