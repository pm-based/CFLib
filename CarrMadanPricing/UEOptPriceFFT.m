function optionPrice = UEOptPriceFFT(spotPrice, strike, rate, TTM, ...
    putFlag, priceModel, modelParams)

% Discretization parameter
Npow=20; N=2^Npow; A=1200;

eta=A/N; v=[0:eta:A*(N-1)/N]; v(1)=1e-22;

lambda=2*pi/(N*eta);              % lambda = 2*pi/A
k=-lambda*N/2+lambda*(0:N-1);     % k = -pi/A*N + lambda*(0, N-1)

% FFT of z_k
CharFunc = @(v) exp(TTM*CharExp(v,priceModel,modelParams));
Z_k = exp(1i*rate*v*TTM).*(CharFunc(v-1i)-1)./(1i*v.*(1i*v+1));

disp('RiskNeutral Check')
CharFunc(-1i)

% Option Price
w = ones(1,N); w(1) = 0.5; w(end) = 0.5;
x = w.*eta.*Z_k.*exp(1i*pi*(0:N-1));
z_k = real(fft(x)/pi);
C = spotPrice*(z_k + max(1-exp(k-rate*TTM),0));
K = spotPrice*exp(k);

% Output
index = find( K>0.1*spotPrice & K < 3*spotPrice );
C = C(index); K = K(index);
plot(K,C);
title( 'Option Price' );
xlabel('Strike');
callPrice = interp1(K,C,strike,'spline');

if putFlag
    optionPrice = callPrice - spotPrice + strike * exp(-rate * TTM);
else
    optionPrice = callPrice;
end


end

function V=CharExp(v,priceModel,modelParams)
switch priceModel
    case 'Merton'
        V = @(u) -(modelParams.sigmaD * u).^2/2 ...
                 + modelParams.lambda ...
                 .* (exp(-(modelParams.sigmaJ*u).^2/2 + 1i*u*modelParams.muJ) -1);

    case 'Kou'
        V = @(u) -(modelParams.sigmaD * u).^2/2 ...
            + 1i*u*modelParams.lambda ...
            .* ( modelParams.p./(modelParams.lambdaP-1i*u) ...
            - (1-modelParams.p)./(modelParams.lambdaN+1i*u) );
end
drift_RN = -V(-1i);
V = drift_RN * 1i*v + V(v);

end