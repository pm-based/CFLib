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

disp('MC price Eu option, Merton model')
[price,CIprice] = UEOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, mertonParams, nSims)
disp('MC price Eu option, Merton model using Antitetich Variable')
[price,CIprice] = UEOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, mertonParams, nSims, true)

%% Kou
priceModel = 'Kou';
kouParams.sigmaD    = 0.7;
kouParams.lambda    = 10;
kouParams.lambdaP   = 15;
kouParams.lambdaN   = 25;
kouParams.p         = 0.6;

disp('MC price Eu option, Kou model')
[price,CIprice] = UEOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, kouParams, nSims)
disp('MC price Eu option, Kou model using Antitetich Variable')
[price,CIprice] = UEOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, kouParams, nSims, true)