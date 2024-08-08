function optionPrice = UEOptBSPricePDE(spotPrice, strike, rate, ...
    sigma, TTM, putFlag, method, PDEParams, barrierType, barrier)

if nargin < 9      % If barrier type is not defined
    barrier = nan;
    barrierType = 'none';
    nMonitoring = 1;
    barrierCheck = @(S_t,barrier) 1;

elseif nargin < 10 % If barrier value is not defined
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

% Convert the putFlag into a numerical value that will be used in the payoff calculation
if putFlag
    putFlag = -1;  % For a put option, payoff involves max(K - S, 0)
else
    putFlag = 1;   % For a call option, payoff involves max(S - K, 0)
end

if strcmp(method, 'Implicit')
    PDEParams.theta = 1;
end
if strcmp(method, 'Explicit')
    PDEParams.theta = 0;
end

% Definitions
drift = (rate - sigma^2/2);

% Truncation
sigma_multiplier = 6;
x_min = drift * TTM - sigma_multiplier * sigma * sqrt(TTM);
x_max = drift * TTM + sigma_multiplier * sigma * sqrt(TTM);

% Discretization
x = linspace(x_min, x_max, PDEParams.N+1)';

dx = (x_max - x_min)/PDEParams.N;
dt = TTM/PDEParams.M;

% Construction of the matrix
comm_part = sigma^2/(dx)^2;

A = PDEParams.theta * (-drift/(2*dx) + comm_part/2);
B = -1/dt - PDEParams.theta * (comm_part + rate);
C = PDEParams.theta * (drift/(2*dx) + comm_part/2);

Matrix = spalloc(PDEParams.N+1, PDEParams.N+1, 3*(PDEParams.N-1)+2);
Matrix(1,1) = 1;
Matrix(PDEParams.N+1, PDEParams.N+1) = 1;

Ah = (PDEParams.theta-1) * (-drift/(2*dx) + comm_part/2);
Bh = -1/dt - (PDEParams.theta-1) * (comm_part + rate);
Ch = (PDEParams.theta-1) * (drift/(2*dx) + comm_part/2);
Matrixh = spalloc(PDEParams.N+1, PDEParams.N+1, 3*(PDEParams.N-1));

for i = 1:PDEParams.N-1
        Matrix(i+1, i:i+2) = [A, B, C];
        Matrixh(i+1, i:i+2) = [Ah, Bh, Ch];
end

% Backward-In-Time Procedure
V = max(putFlag * (spotPrice*exp(x)-strike), 0);
BC = zeros(size(V));

for i = PDEParams.M-1:-1:0
    if putFlag == -1
        BC(1) = strike*exp(-rate*(TTM-i*dt)) - spotPrice*exp(x_min);
    else
        BC(end) = spotPrice*exp(x_max) - strike*exp(-rate*(TTM-i*dt));
    end


    rhs = Matrixh*V + BC;
    V = Matrix \ rhs;
end
figure

plot(spotPrice*exp(x),V);
xlabel('Spot price'); title('Call Price');
optionPrice = interp1(x,V,0,'spline');

end