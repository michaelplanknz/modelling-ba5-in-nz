clear all
close all

dateResults = datenum('07JUL2022');            % label for data files - data up to this date was used for fitting
dateData = datenum('21SEP2022');                % data up to this date will be plotted alongside model output

fLbl = 'scen2';
fIn = sprintf('results/results_%s_%s.mat', fLbl, datestr(dateResults, 'dd-mmm-yyyy') );
load(fIn); 


fNameData = "data/epidataPartial_" + datestr(dateData, 'dd-mmm-yyyy') + ".csv";

figsToSave = [];        % specify a list of figure numbers to save to file
                        

nFilter = length(epiVarsCompact);

admissionsLag = 40;
deathsLag = 10;

% Popn size time series (for 1st realisation only)
popnSize = epiVarsCompact(1).N;


% Get variables for plotting
[newDailyCases0, newDailyCases1, newDailyCases2, newDailyCases3, newDailyCasesr, newDailyHosp0, newDailyHosp1, newDailyHosp2, newDailyHosp3, newDailyHospr, newDailyDeaths0, newDailyDeaths1, newDailyDeaths2, newDailyDeaths3, newDailyDeathsr, Hocc, S, Sw, E1, E2, I, A] = getVarsToPlot(t, epiVarsCompact, parBase);

dataComb = readtable(fNameData);

hospData = readtable("not_for_release/data/MOH_hosp_occupancy.xlsx");

tData = datenum(dataComb.date);
ind = tData <= dateData;
tData = tData(ind);
dataComb = dataComb(ind, :);
hospOccData =  cumsum(dataComb.dailyHosp-dataComb.dailyDischarges);


% Sum over immunity status
newDailyCases = newDailyCases0 + newDailyCases1 + newDailyCases2 + newDailyCases3 + newDailyCasesr;
newDailyHosp = newDailyHosp0 + newDailyHosp1 + newDailyHosp2 + newDailyHosp3 + newDailyHospr;
newDailyDeaths = newDailyDeaths0 + newDailyDeaths1 + newDailyDeaths2 + newDailyDeaths3 + newDailyDeathsr;

% Sum over age bands
newDailyCases0_all = squeeze(sum(newDailyCases0, 2))';
newDailyCases1_all = squeeze(sum(newDailyCases1, 2))';
newDailyCases2_all = squeeze(sum(newDailyCases2, 2))';
newDailyCases3_all = squeeze(sum(newDailyCases3, 2))';
newDailyCasesr_all = squeeze(sum(newDailyCasesr, 2))';
newDailyHosp0_all = squeeze(sum(newDailyHosp0, 2))';
newDailyHosp1_all = squeeze(sum(newDailyHosp1, 2))';
newDailyHosp2_all = squeeze(sum(newDailyHosp2, 2))';
newDailyHosp3_all = squeeze(sum(newDailyHosp3, 2))';
newDailyHospr_all = squeeze(sum(newDailyHospr, 2))';
newDailyDeaths0_all = squeeze(sum(newDailyDeaths0, 2))';
newDailyDeaths1_all = squeeze(sum(newDailyDeaths1, 2))';
newDailyDeaths2_all = squeeze(sum(newDailyDeaths2, 2))';
newDailyDeaths3_all = squeeze(sum(newDailyDeaths3, 2))';
newDailyDeathsr_all = squeeze(sum(newDailyDeathsr, 2))';



% Totals across all age groups & immunity status:
incidenceRel_all = squeeze(1/parBase.tE*sum(E1+E2, 2)./sum(popnSize, 2))';       % incidence relative to popn size 
newDailyCases_all = squeeze(sum(newDailyCases, 2))';
newDailyHosp_all = squeeze(sum(newDailyHosp, 2))';
hospOcc_all = squeeze(sum( Hocc, 2))';
hospCum_all = squeeze(sum( cumsum(newDailyHosp), 2))';
newDailyDeaths_all = squeeze(sum(newDailyDeaths, 2))'; 
CARt = cumsum(newDailyCases)./cumsum(1/parBase.tE*(E1+E2));
CARt_all = cumsum(newDailyCases_all, 2)./cumsum(squeeze(1/parBase.tE*sum(E1+E2, 2))', 2);




% Adjusted version of daily variables for unascertained first infections
newFirstCases = newDailyCases0 + newDailyCases1 + newDailyCases2 + newDailyCases3;
newFirstHosp = newDailyHosp0 + newDailyHosp1 + newDailyHosp2 + newDailyHosp3;
newFirstDeaths = newDailyDeaths0 + newDailyDeaths1 + newDailyDeaths2 + newDailyDeaths3;
newDailyCasesAscert0 = newDailyCases0 .* (1 + (1-CARt).*newDailyCasesr./newFirstCases);
newDailyCasesAscert1 = newDailyCases1 .* (1 + (1-CARt).*newDailyCasesr./newFirstCases);
newDailyCasesAscert2 = newDailyCases2 .* (1 + (1-CARt).*newDailyCasesr./newFirstCases);
newDailyCasesAscert3 = newDailyCases3 .* (1 + (1-CARt).*newDailyCasesr./newFirstCases);
newDailyCasesAscertr = newDailyCasesr .* CARt;
newDailyHospAscert0 = newDailyHosp0 .* (1 + (1-CARt).*newDailyHospr./newFirstHosp);
newDailyHospAscert1 = newDailyHosp1 .* (1 + (1-CARt).*newDailyHospr./newFirstHosp);
newDailyHospAscert2 = newDailyHosp2 .* (1 + (1-CARt).*newDailyHospr./newFirstHosp);
newDailyHospAscert3 = newDailyHosp3 .* (1 + (1-CARt).*newDailyHospr./newFirstHosp);
newDailyHospAscertr = newDailyHospr .* CARt;
newDailyDeathsAscert0 = newDailyDeaths0 .* (1 + (1-CARt).*newDailyDeathsr./newFirstDeaths);
newDailyDeathsAscert1 = newDailyDeaths1 .* (1 + (1-CARt).*newDailyDeathsr./newFirstDeaths);
newDailyDeathsAscert2 = newDailyDeaths2 .* (1 + (1-CARt).*newDailyDeathsr./newFirstDeaths);
newDailyDeathsAscert3 = newDailyDeaths3 .* (1 + (1-CARt).*newDailyDeathsr./newFirstDeaths);
newDailyDeathsAscertr = newDailyDeathsr .* CARt;


% Proportion of new infections that are re-infections
pReinf = squeeze(sum(E2, 2)./sum(E1+E2, 2))';
pReinfAscert = squeeze(sum(newDailyCasesAscertr, 2))'./squeeze(sum(newDailyCases, 2))';

avgImmNPI_all = squeeze(sum(cat(3, epiVarsCompact.avgImmNPI).*(parBase.popDist'), 2))';
avgImmPI_all = squeeze(sum(cat(3, epiVarsCompact.avgImmPI).*(parBase.popDist'), 2))';





% Waning curves
% Generate draws from prior so can plot prior and posterior curves
nParsToFit = 11;
z = rand(nFilter, nParsToFit); 
ThetaPrior = array2table(z, 'VariableNames', {'dateSeed', 'Cstart', 'Cramp', 'rampDays', 'rampStart', 'pTest', 'IFR', 'IHR', 'relaxAlpha', 'MRampDays', 'waneRate'} );
for iRep = 1:nFilter
    parPrior = getPar(ThetaPrior(iRep, :), parBase);
    par = catstruct(parBase, parPrior);
    [tImmCurves, avgImm_2_vs_inf_prior(iRep, :), avgImm_3_vs_inf_prior(iRep, :), avgImm_inf_vs_inf_prior(iRep, :), avgImm_inf_plus2_vs_inf_prior(iRep, :), avgImm_inf_plus3_vs_inf_prior(iRep, :), avgImm_2_vs_sev_prior(iRep, :), avgImm_3_vs_sev_prior(iRep, :), avgImm_inf_vs_sev_prior(iRep, :), avgImm_inf_plus2_vs_sev_prior(iRep, :), avgImm_inf_plus3_vs_sev_prior(iRep, :), avgImm_2_vs_inf_VOC_prior(iRep, :), avgImm_3_vs_inf_VOC_prior(iRep, :), avgImm_inf_vs_inf_VOC_prior(iRep, :), avgImm_inf_plus2_vs_inf_VOC_prior(iRep, :), avgImm_inf_plus3_vs_inf_VOC_prior(iRep, :), avgImm_2_vs_sev_VOC_prior(iRep, :), avgImm_3_vs_sev_VOC_prior(iRep, :), avgImm_inf_vs_sev_VOC_prior(iRep, :), avgImm_inf_plus2_vs_sev_VOC_prior(iRep, :), avgImm_inf_plus3_vs_sev_VOC_prior(iRep, :)] = calcWaneCurves(par); 
    par = catstruct(parBase, parInd(iRep));
    [tImmCurves, avgImm_2_vs_inf(iRep, :), avgImm_3_vs_inf(iRep, :), avgImm_inf_vs_inf(iRep, :), avgImm_inf_plus2_vs_inf(iRep, :), avgImm_inf_plus3_vs_inf(iRep, :), avgImm_2_vs_sev(iRep, :), avgImm_3_vs_sev(iRep, :), avgImm_inf_vs_sev(iRep, :), avgImm_inf_plus2_vs_sev(iRep, :), avgImm_inf_plus3_vs_sev(iRep, :), avgImm_2_vs_inf_VOC(iRep, :), avgImm_3_vs_inf_VOC(iRep, :), avgImm_inf_vs_inf_VOC(iRep, :), avgImm_inf_plus2_vs_inf_VOC(iRep, :), avgImm_inf_plus3_vs_inf_VOC(iRep, :), avgImm_2_vs_sev_VOC(iRep, :), avgImm_3_vs_sev_VOC(iRep, :), avgImm_inf_vs_sev_VOC(iRep, :), avgImm_inf_plus2_vs_sev_VOC(iRep, :), avgImm_inf_plus3_vs_sev_VOC(iRep, :)] = calcWaneCurves(par); 
end




%%
qt = [0.05 0.25 0.75 0.95]; clr = [0 0 1]; ls = '-'; mkr = 'none';
qt90 = [0.05 0.95];
qtiq = [0.25 0.75]; 
greyCol = [0.5 0.5 0.5];
colOrd = colororder;
colOrd = [colOrd; [0 0 0]];
nTraj = 0*min(5, nFilter);
t = datetime(t, 'ConvertFrom', 'datenum');
tData = datetime(tData, 'ConvertFrom', 'datenum');

tPlotRange = [datetime('22JAN2022', 'InputFormat', 'ddMMMyyyy'), ...
    datetime('01OCT2022', 'InputFormat', 'ddMMMyyyy')];









hf = figure(1);
hf.Position = [0 0 1800 800];
subplot(2, 4, 1)
errorShadeFull(t, incidenceRel_all*1e5, qt, clr, ls, mkr )
hold on
if nTraj > 0; plot(t, incidenceRel_all(1:nTraj, :)*1e5, 'Color', greyCol); end
plot(tData, dataComb.dailyInfectionsPer100K, 'o')
xlim(tPlotRange)
grid on
ylabel('new daily infections per 100,000')
title('(a)')
subplot(2, 4, 2)
errorShadeFull(t, newDailyCases_all, qt, clr, ls, mkr )
hold on
if nTraj > 0; plot(t, newDailyCases_all(1:nTraj, :), 'Color', greyCol); end
indFitCases = datenum(tData) <= dateResults-3;        % -3 to allow for fact we are showing 7-day smoothed data
y = smoothdata(dataComb.dailyCases, 'movmean', 7);
plot(tData(indFitCases), y(indFitCases) )
plot(tData(~indFitCases), y(~indFitCases) )
ylm = ylim;
if parInd(1).responseCaseTrig < ylm(2)
   yline(parInd(1).responseCaseTrig, 'k--'); 
end
xlim(tPlotRange)
ylabel('new daily cases')
grid on
title('(b)')
subplot(2, 4, 3)
errorShadeFull(t, newDailyHosp_all, qt, clr, ls, mkr )
hold on
if nTraj > 0; plot(t, newDailyHosp_all(1:nTraj, :), 'Color', greyCol); end
indFitHosp = datenum(tData) <= dateResults-admissionsLag;        
indShowHosp = datenum(tData) > dateResults-admissionsLag & datenum(tData) <= dateData-admissionsLag;
y = smoothdata(dataComb.dailyHosp, 'movmean', 7);
plot(tData(indFitHosp), y(indFitHosp) )
plot(tData(indShowHosp), y(indShowHosp) )
xlim(tPlotRange)
ylabel('new daily admissions for covid')
grid on
title('(c)')
subplot(2, 4, 4)
errorShadeFull(t, newDailyDeaths_all, qt, clr, ls, mkr )
hold on
if nTraj > 0; plot(t, newDailyDeaths_all(1:nTraj, :), 'Color', greyCol); end
indFitDeaths = datenum(tData) <= dateResults-deathsLag;        
indShowDeaths = datenum(tData) > dateResults-deathsLag & datenum(tData) <= dateData-deathsLag;
y = smoothdata(dataComb.dailyDeaths, 'movmean', 7);
plot(tData(indFitDeaths), y(indFitDeaths))
plot(tData(indShowDeaths), y(indShowDeaths))
xlim(tPlotRange)
ylabel('daily deaths')
grid on
title('(d)')
subplot(2, 4, 5)
errorShadeFull(t, cumsum(incidenceRel_all, 2), qt, clr, ls, mkr )
hold on
if nTraj > 0; plot(t, cumsum(incidenceRel_all(1:nTraj, :), 2), 'Color', greyCol); end
y = dataComb.dailyInfectionsPer100K/1e5*7;
ind = ~isnan(y);
plot(tData(ind), cumsum(y(ind)), 'o')
xlim(tPlotRange)
ylabel('cumulative infections per capita')
grid on
title('(e)')
subplot(2, 4, 6)
errorShadeFull(t, cumsum(newDailyCases_all, 2), qt, clr, ls, mkr)
hold on
if nTraj > 0; plot(t, cumsum(newDailyCases_all(1:nTraj, :), 2), 'Color', greyCol); end
y = cumsum(dataComb.dailyCases);
plot(tData(indFitCases), y(indFitCases), '.')
plot(tData(~indFitCases), y(~indFitCases), '.')
xlim(tPlotRange)
ylabel('cumulative cases')
grid on
title('(f)')
subplot(2, 4, 7)
errorShadeFull(t, hospOcc_all , qt, clr, ls, mkr)
hold on
if nTraj > 0; plot(t, hospOcc_all(1:nTraj, :), 'Color', greyCol); end
plot(tData(indFitHosp), hospOccData(indFitHosp), '.')
plot(tData(indShowHosp), hospOccData(indShowHosp), '.')
xlim(tPlotRange)
ylabel('hospital occupancy for covid')
grid on
title('(g)')
subplot(2, 4, 8)
errorShadeFull(t, cumsum(newDailyDeaths_all, 2), qt, clr, ls, mkr )
hold on
if nTraj > 0; plot(t, cumsum(newDailyDeaths_all(1:nTraj, :), 2), 'Color', greyCol); end
y = cumsum(dataComb.dailyDeaths);
plot(dataComb.date(indFitDeaths), y(indFitDeaths), '.')
plot(dataComb.date(indShowDeaths), y(indShowDeaths), '.')
xlim(tPlotRange)
ylabel('cumulative deaths')
grid on
title('(h)')









figure(2);
h = gcf; h.Position = [   680   664   787   314];
subplot(1, 2, 1)
errorShadeFull(t, parBase.R0*Ct, qt, clr, ls, mkr)
hold on
if nTraj > 0; plot(t, parBase.R0*Ct(1:nTraj, :), 'Color', greyCol); end
ylabel('R_{EI}(t)')
xlim(tPlotRange)
ylim([0 5])
grid on
title('(a)')
subplot(1, 2, 2)
errorShadeFull(t, Rt, qt, clr, ls, mkr)
hold on
if nTraj > 0; plot(t, Rt(1:nTraj, :), 'Color', greyCol); end
ylabel('R_{eff}(t)')
yline(1, 'k--');
xlim(tPlotRange)
ylim([0 1.8])
grid on
title('(b)')





hf = figure(3);
hf.Position = [   680   644   733   334];
errorShadeFull(t, pReinf, qt, clr, ls, mkr )
hold on
errorShadeFull(t, pReinfAscert, qt, [1 0 0], ls, mkr )
plot(tData(indFitCases), dataComb.pReinf(indFitCases), '-')
plot(tData(~indFitCases), dataComb.pReinf(~indFitCases), '-', 'Color', [0 0.7 0])
ylabel('proportion of cases that are reinfections')
xlim(tPlotRange)
grid on










figure(4);
subplot(2, 2, 1)
errorShadeFull(t, squeeze(sum(S, 2))', qt, clr, ls, mkr )
hold on
if nTraj > 0; plot(t, squeeze(sum(S(:, :, 1:nTraj), 2))', 'Color', greyCol ); end
ylabel('susceptible')
xlim(tPlotRange)
grid on
subplot(2, 2, 2)
errorShadeFull(t, (squeeze(sum(Sw, 2))./sum(popnSize, 2))', qt, clr, ls, mkr )
hold on
if nTraj > 0; plot(t, (squeeze(sum(Sw(:, :, 1:nTraj), 2))./sum(popnSize, 2))', 'Color', greyCol ); end
ylabel('weighted susceptibility')
xlim(tPlotRange)
grid on
subplot(2, 2, 3)
errorShadeFull(t, squeeze(sum(I+A, 2))', qt, clr, ls, mkr )
hold on
if nTraj > 0; plot(t, squeeze(sum(I(:, :, 1:nTraj)+A(:, :, 1:nTraj), 2))', 'Color', greyCol ); end

ylabel('infectious')
xlim(tPlotRange)
grid on
subplot(2, 2, 4)
errorShadeFull(t, squeeze(1/parBase.tE*sum(E1+E2, 2))', qt, clr, ls, mkr )
hold on
if nTraj > 0; plot(t, squeeze(1/parBase.tE*sum(E1(:, :, 1:nTraj)+E2(:, :, 1:nTraj), 2))', 'Color', greyCol ); end
ylabel('new daily infections')
xlim(tPlotRange)
grid on





figure(5);
set(gcf, 'DefaultAxesColorOrder', colOrd, 'DefaultAxesLineStyleOrder',{'-','--','-.'});
subplot(2, 2, 1)
plot(t, S(:, :, 1))
ylabel('susceptible')
legend('0-4', '5-9', '10-14', '15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49', '50-54', '55-59', '60-64', '65-69', '70-74', '75+')
xlim(tPlotRange)
grid on

subplot(2, 2, 2)
plot(t, I(:, :, 1)+A(:, :, 1))
ylabel('infectious')
xlim(tPlotRange)
grid on

subplot(2, 2, 3)
plot(t, 1/parBase.tE*(E1(:, :, 1)+E2(:, :, 1)))
ylabel('incidence (per day)')
xlim(tPlotRange)
grid on

subplot(2, 2, 4)
plot(t, Sw(:, :, 1)./popnSize)
ylabel('aggregate susceptibility')
grid on




figure(6);
h = gcf; h.Position = [   680         616        1179         362];
h = violinplot(Theta);
ylabel('z_i')
h = gca;
h.XTickLabel = {'date seed', 'R_{EI,1}', 'R_{EI,2}', 'P1-P2 ramp', 'P1 end', 'p_{test,clin}', 'alpha_{IFR}', 'alpha_{IHR}', 'alpha_M', 'M ramp', 'r_w'};


% Prior immunity curves
hf = figure(7);
hf.Position = [     680         342        1119         636];
subplot(2, 2, 1)
errorShadeFull(tImmCurves, avgImm_2_vs_inf_prior, qt90, [0 0 1], ls, mkr)
hold on
errorShadeFull(tImmCurves, avgImm_3_vs_inf_prior, qt90, [1 0 0], ls, mkr)
errorShadeFull(tImmCurves, avgImm_inf_vs_inf_prior, qt90, [0 1 0], ls, mkr)
errorShadeFull(tImmCurves, avgImm_inf_plus2_vs_inf_prior, qt90, [1 0 1], ls, mkr)
errorShadeFull(tImmCurves, avgImm_inf_plus3_vs_inf_prior, qt90, [0 1 1], ls, mkr)
ylabel('immunity to infection with BA.2')
xlabel('time since last immunising event (days)')
grid on
title('(a)')

subplot(2, 2, 2)
errorShadeFull(tImmCurves, avgImm_2_vs_sev_prior, qt90, [0 0 1], ls, mkr)
hold on
errorShadeFull(tImmCurves, avgImm_3_vs_sev_prior, qt90, [1 0 0], ls, mkr)
errorShadeFull(tImmCurves, avgImm_inf_vs_sev_prior, qt90, [0 1 0], ls, mkr)
errorShadeFull(tImmCurves, avgImm_inf_plus2_vs_sev_prior, qt90, [1 0 1], ls, mkr)
errorShadeFull(tImmCurves, avgImm_inf_plus3_vs_sev_prior, qt90, [0 1 1], ls, mkr)
xlabel('time since last immunising event (days)')
ylabel('immunity to severe disease from BA.2')
ylim([0 1])
legend('2 dose', '3 dose', 'prior BA.2', '2 dose + prior BA.2', '3 dose + prior BA.2', 'location', 'SouthWest')
grid on
title('(b)')

subplot(2, 2, 3)
errorShadeFull(tImmCurves, avgImm_2_vs_inf_VOC_prior, qt90, [0 0 1], ls, mkr)
hold on
errorShadeFull(tImmCurves, avgImm_3_vs_inf_VOC_prior, qt90, [1 0 0], ls, mkr)
errorShadeFull(tImmCurves, avgImm_inf_vs_inf_VOC_prior, qt90, [0 1 0], ls, mkr)
errorShadeFull(tImmCurves, avgImm_inf_plus2_vs_inf_VOC_prior, qt90, [1 0 1], ls, mkr)
errorShadeFull(tImmCurves, avgImm_inf_plus3_vs_inf_VOC_prior, qt90, [0 1 1], ls, mkr)
ylabel('immunity to infection with BA.5')
xlabel('time since last immunising event (days)')
grid on
title('(c)')

subplot(2, 2, 4)
errorShadeFull(tImmCurves, avgImm_2_vs_sev_VOC_prior, qt90, [0 0 1], ls, mkr)
hold on
errorShadeFull(tImmCurves, avgImm_3_vs_sev_VOC_prior, qt90, [1 0 0], ls, mkr)
errorShadeFull(tImmCurves, avgImm_inf_vs_sev_VOC_prior, qt90, [0 1 0], ls, mkr)
errorShadeFull(tImmCurves, avgImm_inf_plus2_vs_sev_VOC_prior, qt90, [1 0 1], ls, mkr)
errorShadeFull(tImmCurves, avgImm_inf_plus3_vs_sev_VOC_prior, qt90, [0 1 1], ls, mkr)
xlabel('time since last immunising event (days)')
ylabel('immunity to severe disease from BA.5')
ylim([0 1])
grid on
title('(d)')












figure(9);
plot(t, sum(parBase.doses1, 2)./sum(popnSize, 2), t, sum(parBase.doses2, 2)./sum(popnSize, 2), t, sum(parBase.doses3, 2)./sum(popnSize, 2), t, sum(parBase.doses4plus, 2)./sum(popnSize, 2) )
legend('1st dose', '2nd dose', '3rd dose', '4th dose')
xline(datetime(2022, 7, 11), 'k--', 'HandleVisibility', 'off');
xlim(tPlotRange)
ylim([0 1])
ylabel('proportion vaccniated')
grid on





for iFig = 1:length(figsToSave)
    hf = figure(figsToSave(iFig));
    fOut = append('plots/', fLbl, '_Figure', string(figsToSave(iFig)), '_', datestr(dateData, 'YYYY-mm-DD'), '.png');
    saveas(hf, fOut);
    pause(0.2);     % Sometimes needed to prevent file system problems
end












