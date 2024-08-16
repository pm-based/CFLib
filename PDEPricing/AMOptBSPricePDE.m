function optionPrice = AMOptBSPricePDE(spotPrice, strike, rate, ...
    sigma, TTM, putFlag, PDEParams)

% Convert the putFlag into a numerical value that will be used in the payoff calculation
if putFlag
    putFlag = -1;  % For a put option, payoff involves max(K - S, 0)
else
    putFlag = 1;   % For a call option, payoff involves max(S - K, 0)
end


% Definitions
drift = (rate - sigma^2/2);

% Truncation
sigma_multiplier = 6;
S_min = 0;%spotPrice * exp(drift * TTM - sigma_multiplier * sigma * sqrt(TTM));
S_max = 5;%spotPrice * exp(drift * TTM + sigma_multiplier * sigma * sqrt(TTM));


% Discretization
S = linspace(S_min, S_max, PDEParams.N+1)';

dS = (S_max - S_min)/PDEParams.N;
dt = TTM/PDEParams.M;

% Construction of the matrix
nodes=S(2:end-1);
comm_part = sigma^2/(dS)^2;

A = -rate*nodes / (2*dS) + 0.5*comm_part * nodes.^2;
B = -1/dt - comm_part * nodes.^2 - rate;
C = rate*nodes / (2*dS) + 0.5*comm_part * nodes.^2;

Matrix = spalloc(PDEParams.N+1, PDEParams.N+1, 3*(PDEParams.N-1)+2);
Matrix(1,1) = 1;
Matrix(PDEParams.N+1, PDEParams.N+1) = 1;


for i = 2:PDEParams.N
        Matrix(i, i-1:i+1) = [A(i-1), B(i-1), C(i-1)];
end

% Backward-In-Time Procedure
payoff = max(putFlag * (S-strike), 0);
V = payoff;
BC = zeros(size(V));

for i = PDEParams.M-1:-1:0
    if putFlag == -1
        BC(1) = strike - S_min; % Different wrt EU Check: also different from theory
    else
        BC(end) = S_max - strike;
    end

    rhs = -[0;V(2:end-1);0]/dt+BC;
    V = PSOR(V,Matrix,rhs,payoff);
end
figure

plot(S,V);
xlabel('Spot price'); title('Call Price');
optionPrice = interp1(S,V,spotPrice,'spline');

end

function ynew = PSOR(yold, C, b, p)
maxiter =500; tol = 1e-6; omega = 1.5;
N = length(C); 
ynew = zeros(size(yold));

for j = 1:maxiter
    for i = 1:N
        % exploiting tridiagonal form of the matrix
        if i==1
            z=(b(i)-C(i,i+1)*yold(i+1))/C(i,i);
        elseif i==N
            z=(b(i)-C(i,i-1)*ynew(i-1))/C(i,i);
        else
            z=(b(i)-C(i,i-1)*ynew(i-1)-C(i,i+1)*yold(i+1))/C(i,i);
        end

        ynew(i)=max (p(i), yold(i)+omega*(z-yold(i)));
    end

    if norm( yold-ynew,'inf')<tol
        break
    else
        yold=ynew;
    end

end
end


