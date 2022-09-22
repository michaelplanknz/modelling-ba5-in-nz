function [Rt, gr, Ct] = getGrowthMetrics(t, epiVarsCompact, parBase, parInd)

% Calculate time-varying reproduction number Rt, epidemic growth rate gr, and value of
% control function Ct for an ensemble of model solutions contsined in the
% structure array epiVarsCompasct and with accompanying parameters
% contained in the structure array parInd

nFilter = length(epiVarsCompact);
nDays = length(t);

Rt = zeros(nFilter, nDays);
Ct = zeros(nFilter, nDays);
gr = zeros(nFilter, nDays);
a = 0:0.01:30;
gp1 = exppdf(a, parBase.tE);
gp2 = exppdf(a, parBase.tI);
gp = conv(gp1, gp2);
gp = gp(1:length(a))/trapz(a, gp(1:length(a)));
for iRep = 1:nFilter
   parTemp = catstruct(parBase, parInd(iRep));
   [Rt(iRep, :), Ct(iRep, :)] = calcRt(t, epiVarsCompact(iRep), parTemp);
   gr(iRep, :) = calcGRfromReff(Rt(iRep, :), a, gp);
end
