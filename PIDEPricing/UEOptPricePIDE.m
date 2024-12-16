function optionPrice = UEOptPricePIDE(spotPrice, strike, rate, ...
    TTM, putFlag, priceModel, modelParams, PDEParams, CompleteIntegralFlag, barrierType, barrier)
%UEOPTPRICEPIDE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 9
    CompleteIntegralFlag = false;
elseif nargin < 10      % If barrier type is not defined
    barrier = nan;
    barrierType = 'none';
    nMonitoring = 1;
    barrierCheck = @(S_t,barrier) 1;

elseif nargin < 11 % If barrier value is not defined
    barrier = nan;
    barrierType = 'none';
    nMonitoring = 1;
    barrierCheck = @(S_t,barrier) 1;
    disp('ERROR: no barrier or nMonitoring value specified, no barrier used.')

elseif barrierType == 'DO'
    barrierCheck = @(S_t,barrier) (min(S_t,[],2) > barrier);
    
elseif barrierType == 'DI'
    barrierCheck = @(S_t,barrier) (min(S_t,[],2) < barrier);

elseif barrierType == 'UO'
    barrierCheck = @(S_t,barrier) (max(S_t,[],2) > barrier);

elseif barrierType == 'UI'
    barrierCheck = @(S_t,barrier) (max(S_t,[],2) < barrier);

else
    barrier = nan;
    barrierType = 'none';
    barrierCheck = @(S_t,barrier) 1;
    disp('ERROR: barrier type not allowed, no barrier used.')
end

if strcmpi(priceModel,'NIG') || strcmpi(priceModel, 'VG')
    FiniteActivity = false;
    if ~isfield(modelParams, 'sigma_ext')
        modelParams.sigma_ext = 0;
    end
else
    FiniteActivity = true;
end

% Convert the putFlag into a numerical value that will be used in the payoff calculation
if putFlag
    putFlag = -1;  % For a put option, payoff involves max(K - S, 0)
else
    putFlag = 1;   % For a call option, payoff involves max(S - K, 0)
end

% levy density
k = levyDensity(priceModel, modelParams);

% Truncation
S_min = 0.1*spotPrice;
S_max = 4*spotPrice;
x_min = log(S_min/spotPrice);
x_max = log(S_max/spotPrice);
x = linspace(x_min,x_max,PDEParams.N+1)';
dx = (x_max-x_min)/PDEParams.N;
dt = TTM/PDEParams.M;

% Integral terms
if FiniteActivity
    [B_max, B_min, alpha, lambda_numeric] = Integral_FiniteActiviy(x_max, x_min, PDEParams.N, k);
    sigma = modelParams.sigmaD;
else
    y = linspace(-PDEParams.eps, PDEParams.eps, 2*PDEParams.N);
    sigma2_eps = trapz(y, y.^2.*k(y));
    sigma = sqrt(modelParams.sigma_ext^2 + sigma2_eps);
    if CompleteIntegralFlag
        k=@(y) k(y).*(abs(y)>PDEParams.eps); % AR addition
    end
    [B_max, B_min, alpha, lambda_numeric] = Integral_InfiniteActiviy(x_max, x_min, PDEParams.N, k, PDEParams.eps);

end

if CompleteIntegralFlag
    alpha = 0;
    lambda = 0;
else
    lambda = lambda_numeric;
end

comm_part = sigma^2/dx^2;
drift = (rate - sigma^2/2);


A = PDEParams.theta * ( -(drift - alpha)/(2*dx) + comm_part/2);
B = -1/dt + PDEParams.theta*(-comm_part - (rate+lambda));
C = PDEParams.theta * ( (drift - alpha)/(2*dx) + comm_part/2);

Ah = (PDEParams.theta-1) * ( -(drift - alpha)/(2*dx) + comm_part/2);
Bh = -1/dt + (PDEParams.theta-1)* (-comm_part - (rate+lambda));
Ch = (PDEParams.theta-1) * ( (drift - alpha)/(2*dx) + comm_part/2);

Matrix = spalloc(PDEParams.N+1, PDEParams.N+1, 3*(PDEParams.N+1)+2);
Matrix(1,1) = 1;
Matrix(PDEParams.N+1,PDEParams.N+1) = 1;

Matrixh = spalloc(PDEParams.N+1, PDEParams.N+1, 3*(PDEParams.N+1));

for i = 2:PDEParams.N
    Matrix(i,i-1:i+1) = [A, B, C];
    Matrixh(i,i-1:i+1) = [Ah, Bh, Ch];
end

% Backward loop
V = max(0, spotPrice*exp(x) - strike);
BC = zeros(size(V));

for i = PDEParams.M-1 : -1 : 0
    if FiniteActivity
        I = Integral_FiniteActivity_loop(x,k,B_min,B_max,V,spotPrice,...
            strike*exp(-rate*(TTM-(i+1)*dt)),CompleteIntegralFlag);
    else
        I = Integral_InfiniteActiviy_loop(x,k,B_min,B_max,V,spotPrice,...
            strike*exp(-rate*(TTM-(i+1)*dt)), PDEParams.eps, CompleteIntegralFlag);
    end

    if putFlag == -1
        BC(1) = strike*exp(-rate*(TTM-i*dt)) - spotPrice*exp(x_min);
    else
        BC(end) = spotPrice*exp(x_max) - strike*exp(-rate*(TTM-i*dt));
    end
    rhs = Matrixh*V + BC - I;
    V = Matrix\rhs;
end

figure();
plot(spotPrice*exp(x), V);
xlabel('spot price'); title('option price');

optionPrice = interp1(x,V,0,'spline');


end


%% finity activity support functions

function [B_max, B_min, alpha, lambda_numeric] = Integral_FiniteActiviy(x_max, x_min, N, k)
% can be introduced a if whenever alpha and lambda are not used.
tol = 1e-10;

B_min = x_min;
B_max = x_max;

while k(B_min) > tol
    B_min = B_min-1;
end

while k(B_max) > tol
    B_max = B_max+1;
end

Nq = 2*N;
y = linspace(B_min, B_max, Nq+1);
k = k(y);

alpha = trapz(y, (exp(y)-1).*k);
lambda_numeric = trapz(y,k);
end

function I = Integral_FiniteActivity_loop(x,k,B_min,B_max,V,spotPrice,Kdisc,CompleteIntegralFlag)

% Params
I = zeros(length(x), 1);
y = linspace(B_min, B_max, length(x)); dy = y(2) - y(1);
w = ones(size(y))*dy; w(1)=w(1)/2; w(end)=w(1);

if ~CompleteIntegralFlag
    k = k(y);
    % quadrature computation
    for i = 2 : length(x)-1
        yq = x(i) + y;
        I(i) = sum(w.*f_int(yq,x,V,spotPrice,Kdisc).*k);
    end

else
    k = k(y);
    dx = x(2) - x(1);
    % quadrature computation
    for i = 2 : length(x)-1
        yq = x(i) + y;
        dV = (V(i+1)-V(i-1))/(2*dx);
        I(i) = sum(w.*(f_int(yq,x,V,spotPrice,Kdisc) - V(i) - (exp(y)-1)*dV).*k);
    end
end
end

%% Infinite activity support functions
function [B_max, B_min, alpha, lambda_numeric] = Integral_InfiniteActiviy(x_max, x_min, N, k, eps)
% can be introduced a if whenever alpha and lambda are not used.
tol = 1e-10;

B_min = x_min;
B_max = x_max;

while k(B_min) > tol
    B_min = B_min-1;
end

while k(B_max) > tol
    B_max = B_max+1;
end

Nq = 2*N;

% left side
y = linspace(B_min, -eps, Nq);
k1 = k(y);

alpha = trapz(y, (exp(y)-1).*k1);
lambda_numeric = trapz(y,k1);

% left side
y = linspace(eps, B_max, Nq);
k2 = k(y);

alpha = alpha + trapz(y, (exp(y)-1).*k2);
lambda_numeric = lambda_numeric + trapz(y,k2);
end

function I = Integral_InfiniteActiviy_loop(x,k,B_min,B_max,V,spotPrice,Kdisc, eps, CompleteIntegralFlag)

% Params
I = zeros(length(x), 1);

y1 = linspace(B_min, -eps, length(x)); dy = y1(2) - y1(1);
w1 = ones(size(y1))*dy; w1(1)=w1(1)/2; w1(end)=w1(1);
k1 = k(y1);

y2 = linspace(eps, B_max, length(x)); dy = y2(2) - y2(1);
w2 = ones(size(y2))*dy; w2(1)=w2(1)/2; w2(end)=w2(1);
k2 = k(y2);

% quadrature computation
if CompleteIntegralFlag
    dx = x(2) - x(1);

    for i = 2 : length(x)-1
        dV = (V(i+1)-V(i-1))/(2*dx);
        I(i) = sum(w1.*(f_int(x(i)+y1,x,V,spotPrice,Kdisc) - V(i) - (exp(y1)-1)*dV).*k1) + ...
               sum(w2.*(f_int(x(i)+y2,x,V,spotPrice,Kdisc) - V(i) - (exp(y2)-1)*dV).*k2);
    end
else
    for i = 2 : length(x)-1
        I(i) = sum(w1.*f_int(x(i)+y1,x,V,spotPrice,Kdisc).*k1) + ...
               sum(w2.*f_int(x(i)+y2,x,V,spotPrice,Kdisc).*k2);
    end
end
end

%% Other support functions
function f = f_int(y,x,V,spotPrice,Kdisc)

% Da adattare per il caso put....

% params
f = zeros(size(y));

% boundaries conditions
% lower
id = find(y <= x(1));
f(id) = 0;

% middle
id = find( (y > x(1)).*(y < x(end)) );
f(id) = interp1(x,V,y(id));

% upper
id = find( y > x(end) );
f(id) = spotPrice*exp(y(id)) - Kdisc;

end
