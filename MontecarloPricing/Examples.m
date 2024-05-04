clear
clc

%% UE vanilla option pricing
% simulation params
nSims = 1e5;

% option params
spotPrice = 98;
strike = 100;
rate = 0.0001;
TTM = 1;
putFlag = false;

%% Merton
priceModel = 'Merton';
mertonParams.sigmaD = 0.7;
mertonParams.lambda = 10;
mertonParams.muJ    = -0.01;
mertonParams.sigmaJ = 0.4;


[price,CIprice] = UEOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, mertonParams, nSims);

%% Kou
priceModel = 'Kou';
kouParams.sigmaD    = 0.7;
kouParams.lambda    = 10;
kouParams.lambdaP   = 0.3;
kouParams.lambdaN   = 0.2;
kouParams.p         = 0.4;

[price,CIprice] = UEOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, kouParams, nSims);
