clear
clc

%% UE barrier option pricing
% simulation params
nSims = 1e6;

% Market params
spotPrice = 98;
rate = 0.0001;

% option params
strike = 100;
TTM = 1;
putFlag = false;
barrierType = 'DO';
barrier = 95;
nMonitoring = round(12*TTM);

%% Merton
priceModel = 'Merton';
mertonParams.sigmaD = 0.7;
mertonParams.lambda = 10;
mertonParams.muJ    = -0.01;
mertonParams.sigmaJ = 0.4;

disp('MC price Eu option, Merton model')
[price, CIprice] = UEOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, mertonParams, nSims, 'none', barrierType, barrier, nMonitoring)
width_CI_Standard = diff(CIprice)

disp('MC price Eu option, Merton model using Antitetich Variable')
[price, CIprice] = UEOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, mertonParams, nSims, 'AV', barrierType, barrier, nMonitoring)
width_CI_AV = diff(CIprice)

disp('MC price Eu option, Merton model using Control Variable')
[price, CIprice] = UEOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, mertonParams, nSims, 'CV', barrierType, barrier, nMonitoring)
width_CI_CV = diff(CIprice)


%% Kou
priceModel = 'Kou';
kouParams.sigmaD    = 0.7;
kouParams.lambda    = 10;
kouParams.lambdaP   = 15;
kouParams.lambdaN   = 25;
kouParams.p         = 0.6;

disp('MC price Eu option, Kou model')
[price,CIprice] = UEOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, kouParams, nSims, 'none', barrierType, barrier, nMonitoring)
width_CI_Standard = diff(CIprice)


disp('MC price Eu option, Kou model using Antitetich Variable')
[price,CIprice] = UEOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, kouParams, nSims, 'AV', barrierType, barrier, nMonitoring)
width_CI_AV = diff(CIprice)

disp('MC price Eu option, Kou model using Control Variable')
[price, CIprice] = UEOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, kouParams, nSims, 'CV', barrierType, barrier, nMonitoring)
width_CI_CV = diff(CIprice)
