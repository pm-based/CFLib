function [S_t,t_i] = SimAssetPrice(spotPrice, rate, T, nTimeSteps, nProcesses, priceModel, modelParams)
%SIMASSETPRICE Summary of this function goes here
%   Detailed explanation goes here
switch priceModel
    case 'Merton'
        % Simulate asset prices at maturity using the Merton jump-diffusion model
        [X_t, t_i] = MertonProcess(modelParams.sigmaD, modelParams.lambda, ...
            modelParams.muJ, modelParams.sigmaJ, T, nTimeSteps, nProcesses);
    case 'Kou'
        % Simulate asset prices at maturity using the Kou jump-diffusion model
        [X_t, t_i] = KouProcess(modelParams.sigmaD, modelParams.lambda, ...
                            modelParams.lambdaP, modelParams.lambdaN, ...
                            modelParams.p, T, nTimeSteps, nProcesses);
end

S_t = spotPrice * exp(-rate * t_i + X_t);
end

