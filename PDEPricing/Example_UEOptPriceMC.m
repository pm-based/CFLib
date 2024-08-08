clear
clc

%% UE vanilla option pricing with PDE (BS model)
% Market params
spotPrice = 1;
rate = 0.01;
sigma = 0.4;

% common option params
strike = 1;
TTM = 1;

%% Implicit Euler - Call Option
% option params
putFlag = false;

% PDE params
method = 'Implicit';
PDEParams.M = 2000; % time 
PDEParams.N = 1000; % logprice

optionPrice = UEOptBSPricePDE(spotPrice, strike, rate, ...
    sigma, TTM, putFlag, method, PDEParams)

%% Explicit Euler - Call Option
% option params
putFlag = false;

% PDE params
method = 'Explicit';
PDEParams.M = 2000; % time 
PDEParams.N = 400; % logprice

optionPrice = UEOptBSPricePDE(spotPrice, strike, rate, ...
    sigma, TTM, putFlag, method, PDEParams)



%% Theta Method - Call Option
% option params
putFlag = false;

% PDE params
method = 'Theta';
PDEParams.theta = 0.5;
PDEParams.M = 2000; % time 
PDEParams.N = 1000; % logprice

optionPrice = UEOptBSPricePDE(spotPrice, strike, rate, ...
        sigma, TTM, putFlag, method, PDEParams)


%% Implict Euler
% This part uses the function UEOptBSPricePDE_Implicit. It is only a
% didatic funcion to show the SOR algorith, since the implicit euler is
% already implemented in UEOptBSPricePDE and supports also barrier option.

% option params
putFlag = false;

% PDE params
PDEParams.M = 50; % time 
PDEParams.N = 500; % logprice

SORflag = false;
optionPrice = UEOptBSPricePDE_Implicit(spotPrice, strike, rate, ...
        sigma, TTM, putFlag, PDEParams, SORflag)

%%
SORflag = true;
optionPrice = UEOptBSPricePDE_Implicit(spotPrice, strike, rate, ...
        sigma, TTM, putFlag, PDEParams, SORflag)