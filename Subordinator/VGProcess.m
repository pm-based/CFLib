function Y = VGProcess(muVG,sigmaVG,kVG, T, nTimeSteps, nProcesses)
%VGPROCESS Summary of this function goes here
%   Detailed explanation goes here

Y = zeros(nProcesses,nTimeSteps+1);

dt = T/nTimeSteps;

Z = nrand(nProcesses,nTimeSteps);
dS = kVG * icdf('Gamma', rand(nProcesses,nTimeSteps), dt/k, 1);

for j = 1:nTimeSteps+1
    Y(:,j+1) = Y(:,j) + muVG*dS + sigmaVG*sqrt(dS(:,j)).*Z(:,j);
end

