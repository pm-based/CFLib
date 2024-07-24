%% Example Title
% Summary of example objective

%% Simulation Asset Price: Merton Model
% Market parameters
S0=98; r=0.0001;

% Numerical parameters
T=1; Nsim = 1e1; M=100;

% Model parameters
priceModel = 'Merton';
MertonParams.sigmaD = 0.7;
MertonParams.lambda = 10;
MertonParams.muJ = -0.01;
MertonParams.sigmaJ = 0.4;

[S_t,t_i] = SimAssetPrice(S0, r, T, M, Nsim, priceModel, MertonParams);

[check,~,check_CI]=normfit(S_t(:,end)*exp(-r*T)/S0)  

plot(t_i,S_t(1:5,:))

%% Section 2 Title
%  Market parameters
S0=98; r=0.0001;

% Numerical parameters
T=1; Nsim = 1e1; M=100;

% Model parameters
priceModel = 'Kou';
KouParams.sigmaD = 0.7;
KouParams.lambda = 10;
KouParams.p = 0.6;
KouParams.lambdaP = 15;
KouParams.lambdaN = 25;

[S_t,t_i] = SimAssetPrice(S0, r, T, M, Nsim, priceModel, KouParams);

[check,~,check_CI]=normfit(S_t(:,end)*exp(-r*T)/S0)  

plot(t_i,S_t(1:5,:))