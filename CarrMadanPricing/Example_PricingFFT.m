%% UE vanilla option pricing

% Market params
spotPrice = 98;
rate = 0.0001;

% Option params
strike = 100;
TTM = 1;
putFlag = false;

%% Merton
priceModel = 'Merton';
mertonParams.sigmaD = 0.7;
mertonParams.lambda = 10;
mertonParams.muJ    = -0.01;
mertonParams.sigmaJ = 0.4;

UEOptPriceFFT(spotPrice, strike, rate, TTM, putFlag, priceModel, mertonParams)

%% Kou
priceModel = 'Kou';
kouParams.sigmaD    = 0.7;
kouParams.lambda    = 10;
kouParams.lambdaP   = 15;
kouParams.lambdaN   = 25;
kouParams.p         = 0.6;

UEOptPriceFFT(spotPrice, strike, rate, TTM, putFlag, priceModel, kouParams)