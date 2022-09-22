clear all
close all

% scenarios
nScenarios = 3;
scenarios.vocWaneAmount = [0.19 0.39 0.59];

dateLbl = "07-Jul-2022";            % label for data files 

nReps = 50000;             % 50000      % number of reps to do
pFilter = 0.01;            % 0.01   % proportion of best fitting reps to retain
nParsToFit = 11;                    % number of randomised (fitted) parameters
errWeights = [1 1 1 1 1];         % weights on components of error 



fName = "data/epidataPartial_" + dateLbl + ".csv";
fprintf('   Loading epi data %s\n', fName)
dataComb = readtable(fName);     % read in data for fitting
tMaxData = datenum(dataComb.date(find(~isnan(dataComb.dailyCases), 1, 'last' )));

parBase = getBasePar(tMaxData);     % structure of base parameters that do not change from one realisation to the next - run for period with data (tMaxData)

t = parBase.date0:parBase.date0+parBase.tEnd;
IC = getIC(parBase);        % generate an initial condition to find out how many independent variables the ODE has
options = odeset('NonNegative', ones(size(IC))' );


if nReps > 1
   z = rand(nReps, nParsToFit);             % random unifornm values used by getPar to randomise parameter values
   nFilter = ceil(nReps*pFilter);
else
    z = 0.5*ones(1, nParsToFit);        % call getPar for default (median) prior parameter values by setting Theta = 0.5    
    nFilter = 1;
end

% Generate nReps sets of random parameters on the U[0,1] scale (these will
% be transformed onto the prior for the relevant parameter in getPar()
Theta = array2table(z, 'VariableNames', {'dateSeed', 'Cstart', 'Cramp', 'rampDays', 'rampStart', 'pTest', 'IFR', 'IHR', 'relaxAlpha', 'MRampDays', 'waneRate'} );

dist = zeros(nReps, length(errWeights));
fprintf('Running ABC:\n')
parfor iRep = 1:nReps              % parfor
   % parameter structure for ith realisation
   parTemp = getPar(Theta(iRep, :), parBase);
   % temporary structure with merged fields from parBase and parInd
   % because parfor won't allow index variables in an anonymous
   % function:
   parTemp = catstruct(parBase, parTemp);       
    % Initial condition
    IC = getIC(parTemp);
    % solve ODE
    [~, Y] = ode45(@(t, y)myODEs2(t, y, parTemp), t, IC, options);
    % retain a compact structure of variables rather than the full Y
    epiVarsCompactTemp = extractEpiVarsCompact(t, Y, parTemp);
    % calculate error function to data
    dist(iRep, :) = calcError(t, epiVarsCompactTemp, dataComb, parTemp);
end
% Filter reps
[~, keepInd] = mink(sum(errWeights.*dist, 2), nFilter);
Theta = Theta(keepInd, :);
distKeep = dist(keepInd, :);




% Run full length simulations (for each scenario specified) with parameter
% combinations from retained trajectories
parBase = getBasePar(nan);     % structure of base parameters that do not change from one realisation to the next - call with nan to run for full period
t = parBase.date0:parBase.date0+parBase.tEnd;

for iScenario = 1:nScenarios
    fprintf('Scenario %i/%i\n', iScenario, nScenarios)
    parfor iRep = 1:nFilter              % parfor
       % parameter structure for ith realisation
       parInd(iRep) = getPar(Theta(iRep, :), parBase);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Fix specific scenario parameters
       parInd(iRep).vocWaneAmount = scenarios.vocWaneAmount(iScenario);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       % temporary structure with merged fields from parBase and parInd
       % because index variables not allowed in an anonymous
       % function within a parfor loop:
       parTemp = catstruct(parBase, parInd(iRep));       
        % Initial condition
        IC = getIC(parTemp);
        % solve ODE
        [~, Y] = ode45(@(t, y)myODEs2(t, y, parTemp), t, IC, options);
        % retain a compact structure of variables rather than the full Y
        epiVarsCompact(iRep) = extractEpiVarsCompact(t, Y, parTemp);
     end

    % Calculate Rt and growth rate for each accepted realisation
    [Rt, gr, Ct] = getGrowthMetrics(t, epiVarsCompact, parBase, parInd);

    % Save scenario results
    fOut = sprintf('results/results_scen%d_%s.mat', iScenario, dateLbl);
    save(fOut, 't', 'epiVarsCompact', 'Theta', 'dist', 'distKeep', 'Rt', 'Ct', 'gr', 'parBase', 'parInd');

end












