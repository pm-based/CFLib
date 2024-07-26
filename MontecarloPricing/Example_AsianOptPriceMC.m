clear
clc

%% Asian option pricing (floating strike)
% simulation params
nSims = 1e5;

% Market params
spotPrice = 200;
rate = 0.0001;

% Contract params
TTM = 0.5;
nMonitoring = round(52*TTM); 
putFlag = true;

%% Merton
priceModel = 'Merton';
mertonParams.sigmaD = 0.7;
mertonParams.lambda = 10;
mertonParams.muJ    = -0.01;
mertonParams.sigmaJ = 0.4;

disp('MC price Asian option, Merton model')
[price, CIprice] = AsianOptPriceMC(spotPrice, rate, TTM, putFlag, priceModel, mertonParams, nSims, nMonitoring)
width_CI_Standard = diff(CIprice)

disp('MC price Asian option, Merton model using Antitetich Variable')
[price, CIprice] = AsianOptPriceMC(spotPrice, rate, TTM, putFlag, priceModel, mertonParams, nSims, nMonitoring, 'AV')
width_CI_AV = diff(CIprice)

disp('MC price Asian option, Merton model using Control Variable')
[price, CIprice] = AsianOptPriceMC(spotPrice, rate, TTM, putFlag, priceModel, mertonParams, nSims, nMonitoring, 'CV')
width_CI_CV = diff(CIprice)

%% Kou
priceModel = 'Kou';
kouParams.sigmaD    = 0.7;
kouParams.lambda    = 10;
kouParams.lambdaP   = 15;
kouParams.lambdaN   = 25;
kouParams.p         = 0.6;

disp('MC price Asian option, Kou model')
[price,CIprice] = AsianOptPriceMC(spotPrice, rate, TTM, putFlag, priceModel, kouParams, nSims, nMonitoring)
width_CI_Standard = diff(CIprice)

disp('MC price Asian option, Kou model using Antitetich Variable')
[price,CIprice] = AsianOptPriceMC(spotPrice, rate, TTM, putFlag, priceModel, kouParams, nSims, nMonitoring, 'AV')
width_CI_AV = diff(CIprice)

disp('MC price Asian option, Kou model using Control Variable')
[price, CIprice] = AsianOptPriceMC(spotPrice, rate, TTM, putFlag, priceModel, kouParams, nSims, nMonitoring, 'CV')
width_CI_CV = diff(CIprice)
