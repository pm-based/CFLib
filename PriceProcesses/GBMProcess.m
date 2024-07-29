function [X_t, t_i, X_t_AV] = GBMProcess(sigma, T, nTimeSteps, nProcesses, flagAV)
[X_t, t_i, X_t_AV] = MertonProcess(sigma, 0, 0, 0, T, nTimeSteps, nProcesses, flagAV);
end