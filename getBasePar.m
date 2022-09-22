function par = getBasePar(tMax)

% Create a structure of fixed base parameters that do not change from run to run

% Filename label - set to "national" for national model
regionLabel = "national";

% Simulation start date and time period
par.date0 = datenum('05MAR2021');       % 5 MAR 2021 (so that full period of vaccine rollout is simulated to get correct waning dynamics)
if isnan(tMax)
    par.tEnd = 2 * 360;         % no tMax specified run for full period
else
    par.tEnd = tMax-par.date0;  % run for specified period (typically period for which data is available)
end
%------------- SEIR parameters --------------
par.R0 = 3.25;  % Nominal R0 (REI is R0 times Ct, as specified in getPar )
par.cSub = 0.5; % Relative infectiousness of subclinicals
par.tE = 1;     % Latent period
par.tI = 2.3;   % Infectious period



%------------- Specify Population Structure and contact matrix -------------
par.nAgeGroups = 16;

par.nSusComp = 14;
par.nVaxComp = 3;
par.nCaseComp = 3;
par.nHospComp = 5;
par.nDeathComp = 6;

fs = sprintf('data/popsize_%s.xlsx', regionLabel);
par.popCount = readPopnData(fs, par.nAgeGroups);
par.totalPopSize = sum(par.popCount);
par.popDist = par.popCount/sum(par.popCount); 

[par.C_detBal, par.popDistBench] = getC_detBal(par.nAgeGroups);



% ------------------ Population dynamics parameters ---------------------
% Demographic parameters birth, death, ageing - set all of these to zero to
% just have a static population
[Mu, b] = getDemogPars();
par.popnDeathRate = Mu;
par.popnBirthRate = b;
par.popnAgeingRate = 1/(5*365.25);




%------------- Get vaccine data -------------------
par.vaccImmDelay = 14;     % delay from vaccination to immunity
vaxDataDate = datenum("04JUL2022");     % date of most recent vaccine data download
[dates, par.doses1, par.doses2, par.doses3, par.doses4plus] = getVaccineData(regionLabel, vaxDataDate, par);




smoothWindow = 56;      
par.nDoses1Smoothed = [ zeros(1, par.nAgeGroups); diff(smoothdata(par.doses1, 'movmean', smoothWindow));];
par.nDoses2Smoothed = [ zeros(1, par.nAgeGroups); diff(smoothdata(par.doses2, 'movmean', smoothWindow))];
par.nDoses3Smoothed = [ zeros(1, par.nAgeGroups); diff(smoothdata(par.doses3, 'movmean', smoothWindow))];
par.nDoses4Smoothed = [ zeros(1, par.nAgeGroups); diff(smoothdata(par.doses4plus, 'movmean', smoothWindow))];

