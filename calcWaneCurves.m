function [tSpan, avgImm_2_vs_inf, avgImm_3_vs_inf, avgImm_inf_vs_inf, avgImm_inf_plus2_vs_inf, avgImm_inf_plus3_vs_inf, avgImm_2_vs_sev, avgImm_3_vs_sev, avgImm_inf_vs_sev, avgImm_inf_plus2_vs_sev, avgImm_inf_plus3_vs_sev, avgImm_2_vs_inf_VOC, avgImm_3_vs_inf_VOC, avgImm_inf_vs_inf_VOC, avgImm_inf_plus2_vs_inf_VOC, avgImm_inf_plus3_vs_inf_VOC, avgImm_2_vs_sev_VOC, avgImm_3_vs_sev_VOC, avgImm_inf_vs_sev_VOC, avgImm_inf_plus2_vs_sev_VOC, avgImm_inf_plus3_vs_sev_VOC] = calcWaneCurves(par) 

% Calculate waning curves by calculating the fraction of cohort immunised
% (by vaccine or by infection) at time 0 who are in each susceptible
% compartment at time t

wR = par.waneRate_RtoS;
wS = par.waneRate_StoS;

icur = [1; 0; 0; 0];        % Initial condition for immunity curves with no prior infection (have no initial recovered state)
icr = [1; 0; 0; 0; 0];       % Initial condition for immunity curves with prior infection (have an initial recovered state)
dest2 = par.waneNet_RtoS(3, 11:14)';    % Probability of moving to states 11-14 following recovery for 2-dosed people
dest0 = par.waneNet_RtoS(1, 11:14)';    % Probability of moving to states 11-14 following recovery for 0/1-dosed people

tSpan = 0:360;

% solve cohort equations to find out fraction of cohort in each of R, S11,
% S12, S13, S14 as a function of time
myRHS = @(t, y)( [-wS*y(1); wS*(y(1)-y(2)); wS*(y(2)-y(3)); wS*y(3) ] );
[t, Yur] = ode23(myRHS, tSpan, icur);
myRHS = @(t, y)( [-wR*y(1); wR*y(1)-wS*y(2); wS*(y(2)-y(3)); wS*(y(3)-y(4)); wS*y(4) ] );
[t, Y3r] = ode23(myRHS, tSpan, icr);
myRHS = @(t, y)( [-wR*y(1); wR*dest2(1)*y(1)-wS*y(2); wR*dest2(2)*y(1)+wS*(y(2)-y(3)); wR*dest2(3)*y(1)+wS*(y(3)-y(4)); wR*dest2(4)*y(1)+wS*y(4) ] );
[t, Y2r] = ode23(myRHS, tSpan, icr);
myRHS = @(t, y)( [-wR*y(1); wR*dest0(1)*y(1)-wS*y(2); wR*dest0(2)*y(1)+wS*(y(2)-y(3)); wR*dest0(3)*y(1)+wS*(y(3)-y(4)); wR*dest0(4)*y(1)+wS*y(4) ] );
[t, Y0r] = ode23(myRHS, tSpan, icr);

VEi_2 = [par.VEi(3:6)];
VEi_3 = [par.VEi(7:10)];
VEi_r = [1, par.VEi(11:14)];        % previously infection immunity curves have 100% "VE" in the initial recovered state
VEh_2 = [par.VEh(3:6)];
VEh_3 = [par.VEh(7:10)];
VEh_r = [1, par.VEh(11:14)];

% Map cohort fractions onto average immunity
avgImm_2_vs_inf = sum(Yur.*VEi_2, 2);
avgImm_3_vs_inf = sum(Yur.*VEi_3, 2);
avgImm_inf_vs_inf = sum(Y0r.*VEi_r, 2);
avgImm_inf_plus2_vs_inf = sum(Y2r.*VEi_r, 2);
avgImm_inf_plus3_vs_inf = sum(Y3r.*VEi_r, 2);
avgImm_2_vs_sev = sum(Yur.*VEh_2, 2);
avgImm_3_vs_sev = sum(Yur.*VEh_3, 2);
avgImm_inf_vs_sev = sum(Y0r.*VEh_r, 2);
avgImm_inf_plus2_vs_sev = sum(Y2r.*VEh_r, 2);
avgImm_inf_plus3_vs_sev = sum(Y3r.*VEh_r, 2);









% To get immunity curves for BA5, apply an additional "waning period"
% equivalent to time tVOCwane = vocWaneAmount/waneRate_StoS
tVOCwane = par.vocWaneAmount/par.waneRate_StoS;
myRHS = @(t, y)( [-wR*y(1); wR*y(1)-wS*y(2); wS*(y(2)-y(3)); wS*(y(3)-y(4)); wS*y(4) ] );
[~, Y] = ode23(myRHS, [0 tVOCwane], icr);
[t, Y3r] = ode23(myRHS, tSpan, Y(end, :)');
myRHS = @(t, y)( [-wR*y(1); wR*dest2(1)*y(1)-wS*y(2); wR*dest2(2)*y(1)+wS*(y(2)-y(3)); wR*dest2(3)*y(1)+wS*(y(3)-y(4)); wR*dest2(4)*y(1)+wS*y(4) ] );
[~, Y] = ode23(myRHS, [0 tVOCwane], icr);
[t, Y2r] = ode23(myRHS, tSpan, Y(end, :)');
myRHS = @(t, y)( [-wR*y(1); wR*dest0(1)*y(1)-wS*y(2); wR*dest0(2)*y(1)+wS*(y(2)-y(3)); wR*dest0(3)*y(1)+wS*(y(3)-y(4)); wR*dest0(4)*y(1)+wS*y(4) ] );
[~, Y] = ode23(myRHS, [0 tVOCwane], icr);
[t, Y0r] = ode23(myRHS, tSpan, Y(end, :)');


VEi_2 = [par.VEi_VOC(3:6)];
VEi_3 = [par.VEi_VOC(7:10)];
VEi_r = [1, par.VEi_VOC(11:14)];        % previously infection immunity curves have 100% "VE" in the initial recovered state
VEh_2 = [par.VEh_VOC(3:6)];
VEh_3 = [par.VEh_VOC(7:10)];
VEh_r = [1, par.VEh_VOC(11:14)];

avgImm_2_vs_inf_VOC = sum(Yur.*VEi_2, 2);
avgImm_3_vs_inf_VOC = sum(Yur.*VEi_3, 2);
avgImm_inf_vs_inf_VOC = sum(Y0r.*VEi_r, 2);
avgImm_inf_plus2_vs_inf_VOC = sum(Y2r.*VEi_r, 2);
avgImm_inf_plus3_vs_inf_VOC = sum(Y3r.*VEi_r, 2);
avgImm_2_vs_sev_VOC = sum(Yur.*VEh_2, 2);
avgImm_3_vs_sev_VOC = sum(Yur.*VEh_3, 2);
avgImm_inf_vs_sev_VOC = sum(Y0r.*VEh_r, 2);
avgImm_inf_plus2_vs_sev_VOC = sum(Y2r.*VEh_r, 2);
avgImm_inf_plus3_vs_sev_VOC = sum(Y3r.*VEh_r, 2);

