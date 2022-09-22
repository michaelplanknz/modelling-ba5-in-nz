function dist = calcError(t, epiVarsCompact, dataComb, par)

% calculate error function for model output (epiVarsCompact) and data
% (dataComb)

dist  = zeros(1, 5);

minCaseFitDate = datenum('01MAR2022');  % only fit to cases after 1 march
minOtherFitDate = datenum('01FEB2022'); % only fit to other data after 1 feb
deathLag = 10;                      % exclude this number of data points at the end of the time series for deaths to account for date of death to report lag
hospLag = 40;                       % exclude this number of data points at the end of the time series for new admissions to account for report lag
latestData = datenum(dataComb.date( find(~isnan(dataComb.dailyCases), 1, 'last') ));

nDays = length(t);

% error functions will be a sum of squares of log( (x_model + sml)/(x_data + sml))
% sml is needed to avoid divide by 0 errors
% choose sml to be small relative to typical minimum values of each of ther
% 4 fitted variables (daily incidence per capita, daily cases, hospital
% occupancy, daily deaths)
sml_incid = 5e-6;
sml_cases = 10;
sml_hospAdm = 0.1;
sml_deaths = 0.01;
sml_age_frac = 5e-5;

% Get variables for comparing to data
[newDailyCases0, newDailyCases1, newDailyCases2, newDailyCases3, newDailyCasesr, newDailyHosp0, newDailyHosp1, newDailyHosp2, newDailyHosp3, newDailyHospr, newDailyDeaths0, newDailyDeaths1, newDailyDeaths2, newDailyDeaths3, newDailyDeathsr, ~, ~, ~, E1, E2] = getVarsToPlot(t, epiVarsCompact, par);


% Sum over immunity status
newDailyCases = newDailyCases0+newDailyCases1+newDailyCases2+newDailyCases3+newDailyCasesr;
newDailyHosp = newDailyHosp0+newDailyHosp1+newDailyHosp2+newDailyHosp3+newDailyHospr;
newDailyDeaths = newDailyDeaths0+newDailyDeaths1+newDailyDeaths2+newDailyDeaths3+newDailyDeathsr;

% Totals across all age groups:
incidenceRel_all = (1/par.tE * sum(E1+E2, 2)./sum(epiVarsCompact.N, 2));
newDailyCases_all = sum(newDailyCases, 2);
newDailyHosp_all = sum(newDailyHosp, 2);
newDailyDeaths_all = sum(newDailyDeaths, 2); 

% dist(1) - daily cases
ind1 = ismember(t, datenum(dataComb.date)) & t >= minCaseFitDate ;
ind2 = ismember(datenum(dataComb.date), t) & datenum(dataComb.date) >= minCaseFitDate;
assert(isequal(t(ind1)', datenum(dataComb.date(ind2))));
y1 = newDailyCases_all(ind1);
y2 = smoothdata(dataComb.dailyCases, 'movmean', 7);
y2 = y2(ind2);
dist(1) = nanmean( (log(y1+sml_cases) - log(y2+sml_cases)).^2  );



% dist(2) - daily deaths
ind1 = ismember(t, datenum(dataComb.date)) & t >= minOtherFitDate                      & t <= latestData-deathLag;
ind2 = ismember(datenum(dataComb.date), t) & datenum(dataComb.date) >= minOtherFitDate & datenum(dataComb.date) <= latestData-deathLag;
assert(isequal(t(ind1)', datenum(dataComb.date(ind2))));
y1 = newDailyDeaths_all(ind1);
y2 = smoothdata(dataComb.dailyDeaths, 'movmean', 7);
y2 = y2(ind2);
dist(2) = nanmean( (log(y1+sml_deaths) - log(y2+sml_deaths)).^2  );


% dist(3) - daily new admissions
ind1 = ismember(t, datenum(dataComb.date)) & t >= minOtherFitDate                      & t <= latestData-hospLag;
ind2 = ismember(datenum(dataComb.date), t) & datenum(dataComb.date) >= minOtherFitDate & datenum(dataComb.date) <= latestData-hospLag;
assert(isequal(t(ind1)', datenum(dataComb.date(ind2))));
y1 = newDailyHosp_all(ind1);
y2 = smoothdata(dataComb.dailyHosp, 'movmean', 7);
y2 = y2(ind2);
dist(3) = nanmean( (log(y1+sml_hospAdm) - log(y2+sml_hospAdm)).^2  );



% dist(4) - daily incidence per capita
y1 = incidenceRel_all(ind1);
y2 = dataComb.dailyInfectionsPer100K(ind2)/1e5;
assert(isequal(t(ind1)', datenum(dataComb.date(ind2))));
dist(4) = nanmean( (log(y1+sml_incid) - log(y2+sml_incid)).^2  );


% dist(5) - age breakdown of cases (fraction 60+)
ind1 = ismember(t, datenum(dataComb.date)) & t >= minCaseFitDate ;
ind2 = ismember(datenum(dataComb.date), t) & datenum(dataComb.date) >= minCaseFitDate;
assert(isequal(t(ind1)', datenum(dataComb.date(ind2))));
newDailyCases10 = combineBands(newDailyCases);
y1 = newDailyCases10(ind1, :);
% 60+ is the 7th and 8th age bands
y1 = sum(y1(:, 7:8) ,2)./sum(y1, 2);
y2 = dataComb.pOver60(ind2);
dist(5) = nanmean(nanmean( (log(y1+sml_age_frac) - log(y2+sml_age_frac)).^2  ));

