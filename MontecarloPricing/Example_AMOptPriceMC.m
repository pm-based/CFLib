clear
clc

%% American option pricing
% simulation params
nSims = 1e6;

% Market params
spotPrice = 1;
rate = 0.04;

% Contract params
strike = 1;
TTM = 1;
nMonitoring = 40; 
putFlag = true;

%% GBM
priceModel = 'GBM';
BSParams.sigma = 0.4;

disp('MC price AM option, BS model')
[price, CIprice] = AMOptPriceMC(spotPrice, strike, rate, TTM, putFlag, priceModel, BSParams, nSims, nMonitoring)
width_CI_Standard = diff(CIprice)


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
