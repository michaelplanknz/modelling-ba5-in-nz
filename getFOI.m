function [FOI] = getFOI(t, N, I, A, C, VEt, par)

% Calculate force of infection in each susceptible compartment and each age group at time t

dailyCaseRate = 2/par.tLatentToTest * sum(C(:, 2), 1);
Ct = getTimeDepCt(t, dailyCaseRate, par);

% Effect of seed case on FOI
Iseed = par.initialExp * normpdf(t-par.dateSeed, 0, par.seedDur/4);

[~, NGMclin] = getNGMtimeDep(t, par);
FOI = par.R0*Ct/par.tI * (NGMclin * (Iseed*par.tI + sum((1-VEt).*(I + par.cSub*A), 2)) ) ./ N;


