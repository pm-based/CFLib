function levyDensity = levyDensity(priceModel,modelParams)
%LEVYDENSITY Summary of this function goes here
%   Detailed explanation goes here


switch priceModel
    case 'Merton'
        levyDensity = mertonDensity(modelParams.lambda, modelParams.muJ, ...
            modelParams.sigmaJ);

    case 'Kou'
        levyDensity = kouDensity(modelParams.lambda, modelParams.lambdaP, ...
            modelParams.lambdaN, modelParams.p);

    case 'NIG'
        levyDensity = NIGDensity(modelParams.sigmaJ, modelParams.thetaJ, ...
            modelParams.kJ);

    case 'VG'
        levyDensity = VGDensity(modelParams.sigmaJ, modelParams.thetaJ, ...
            modelParams.kJ);
end
end


function mertonDensity = mertonDensity(lambda, muJ, sigmaJ)
mertonDensity = @(y) lambda/(sigmaJ*sqrt(2*pi)) * exp(- (y-muJ).^2/(2*sigmaJ^2));
end

function kouDensity = kouDensity(lambda, lambdaP, lambdaN, p)
kouDensity = @(y) p*lambda*lambdaP * exp(-lambdaP*y) * (y>0) ...
                    + (1-p)*lambda*lambdaP * exp(-lambdaN*abs(y)) * (y<0);
end

function NIGDensity = NIGDensity(sigmaJ, thetaJ, kJ)
A = thetaJ/sigmaJ^2;
N = sqrt(thetaJ^2 + sigmaJ^2/kJ);
B = N/sigmaJ^2;
C = N/(pi*sigmaJ*sqrt(kJ));

NIGDensity = @(y) C./abs(y) .* exp(A*y) .* besselk(1,B*abs(y));
end

function VGDensity = VGDensity(sigmaJ, thetaJ, kJ)
A = thetaJ/sigmaJ^2;
B = sqrt(thetaJ^2 + sigmaJ^2/kJ)/sigmaJ^2;

VGDensity = @(y) 1./(abs(y)*kJ) .* exp(A*y-B*abs(y));
end