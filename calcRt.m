function [Rt, Ct] = calcRt(t, epiVarsCompact, par)

% Calculate the effectvie reproduction number in the model at time t by
% calculating the bext generation matrix (NGM) for a fully susceptible population, the relative
% contact rate Ct, and the aggregate susceptibility at time t

nDays = length(t);
domEig = zeros(1, nDays);


totCumCases = sum(epiVarsCompact.C0+epiVarsCompact.C1+epiVarsCompact.C2+epiVarsCompact.C3+epiVarsCompact.Cr, 2);
dailyCaseRate = [totCumCases(1), diff(totCumCases)'];
Ct = getTimeDepCt(t, dailyCaseRate, par);


for iDay = 1:nDays
    NGM = getNGMtimeDep(t(iDay), par);
    M = NGM.*(epiVarsCompact.Sw(iDay, :)./epiVarsCompact.N(iDay, :))' ;
    domEig(iDay) = eigs(M, 1);
end

Rt = par.R0 * Ct .* domEig;
