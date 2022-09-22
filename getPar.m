function par = getPar(Theta, parBase)


% Functions to add a zero-centred random variable to selected parameters
% according to random deviates [0,1] specified in input parameter Theta
plusMinus = @(z, r)(2*r*z - r);                     % uniform perturbation between +/- r, z is a random deviate [0,1]
plusMinusNorm = @(z, r)( r*norminv(z) );             % normal pertrubation with s.d. r, z is a random deviate [0,1]
plusMinusInt = @(z, r)(floor((2*r+1)*z) - r );      % uniform perturbation on integers between +/- r, z is a random deviate [0,1]


%--------------------- Overall R0 and seed Parameters --------------------
par.dateSeed = datenum('19JAN2022') + plusMinusInt(Theta.dateSeed, 3);
par.seedDur = 7;
par.initialExp = 0.0001*parBase.popCount;


% -------------------- Control function parameters -----------------------
% starting value of Ct
par.Cstart = 0.68 + plusMinus(Theta.Cstart, 0.06);   
 
% value to ramp to after March peak
par.Cramp = 1.2 + plusMinus(Theta.Cramp, 0.3);      
par.rampDays = 55 + plusMinusInt(Theta.rampDays, 20);                               
par.rampStart = datenum('10MAR2022') + plusMinusInt(Theta.rampStart, 5);     

% can be used to add future uncertainty in control function
par.Cramp2 = par.Cramp;
par.rampStart2 = datenum('01SEP2022');
par.rampDays2 = 30;


% dynamic response
par.responseStartDate = datenum('15JUN2022');
par.responseStartWindow = 7;                        % these parameters just ensure the response doesn't apply to the 1st wave
par.responseCt = 0.85*par.Cramp;                     % reduction in contact during control period
par.responseCaseTrig = 1e7;                         % response triggered with a (midpoint) of this number of daily cases
par.responseTrigSmoothWindow = 0.1;                 % change in control function takes place smoothly with s.d. = this number x responseCaseTrig


% Use parameters values to construct control function Ct
t = parBase.date0:max(parBase.date0+parBase.tEnd, par.rampStart2+par.rampDays2+1);     % this might be an overly long time period (for purposes of applying ramp changes in Ct) but will be trunctaed to approproate length below 
par.Ct = par.Cstart*ones(1, parBase.tEnd+1);
ti = find(t == par.rampStart);
par.Ct(ti:ti+par.rampDays-1) = linspace(par.Cstart, par.Cramp, par.rampDays);  
par.Ct(ti+par.rampDays:end) = par.Cramp;
ti = find(t == par.rampStart2);
par.Ct(ti:ti+par.rampDays2-1) = linspace(par.Cramp, par.Cramp2, par.rampDays2);  
par.Ct(ti+par.rampDays2:end) = par.Cramp2;

t = parBase.date0:(parBase.date0+parBase.tEnd);
par.Ct = par.Ct(1:length(t));                   % truncate par.Ct if it slonger than required simulation period length(t)

% --------------------- Testing and lag parameters -----------------------
par.tLatentToTest = 4;              % onset of infectiousness to test
par.tTestToHosp = 1;                % test to hospital admission
par.tLOS = [ 2.0000    2.0000    2.0000    2.0000    2.0000    2.0000    2.6700    3.3400    4.0100    4.6800    5.3500    6.0200    6.6900    7.3600    8.0300    8.7000]';
par.tDeath = 14;                    % Time from admission to death 

par.pTestClin = 0.55 + plusMinus(Theta.pTest, 0.2);            % testing probability clinical
par.pTestSub = 0.4*par.pTestClin;                            % testing probability subclinical


%------------- Disease Rate Data --------------
par.IFRmult = 0.8 * ( 1 + plusMinus(Theta.IFR, 0.5) );  % IFR multiplier
par.IHRmult = 0.3 * ( 1 + plusMinus(Theta.IHR, 0.5));    % IHR multiplier

par.pClin = [0.5440, 0.5550, 0.5770, 0.5985, 0.6195, 0.6395, 0.6585, 0.6770, 0.6950, 0.7117, 0.7272, 0.7418, 0.7552, 0.7680, 0.7800, 0.8008]'; % Fraser group
par.ui = [0.4000, 0.3950, 0.3850, 0.4825, 0.6875, 0.8075, 0.8425, 0.8450, 0.8150, 0.8050, 0.8150, 0.8350, 0.8650, 0.8450, 0.7750, 0.7400 ];    % Davies relative suscewptibility
[IHR0, ~, IFR0] = getHerreraRatesOmi();
par.IFR = par.IFRmult*IFR0;
par.IHR = par.IHRmult*IHR0;


% --------------- contact matrix adjustment parameters ------------------
ageBlockSizes = [3 2 2 3 2 4];      % divide the contact matrix up into blocks - this vector specifies the number of 5-year age classes in each "block"
nAgeBlocks = length(ageBlockSizes);
Cw1 = [1.1 0.7 0.55 0.45 0.45 0.5;  % weighting matrix for initial contact matrix
      0    1.2 0.7  0.5  0.5  0.3;
      0    0   1.1  0.5  0.5  0.5;
      0    0   0    0.15 0.15 0.45;
      0    0   0    0    0.15 0.45;
      0    0   0    0    0    0.15];
par.relaxAlpha = 0.4 + plusMinus(Theta.relaxAlpha, 0.4);                         % Amount by which contract matrix relaxes back to Prem (0=not at all, 1=fully)
Cw2 = (1-par.relaxAlpha)*Cw1 + par.relaxAlpha*triu(ones(nAgeBlocks)) ;
Cw1 = Cw1+triu(Cw1, 1)';                % make weights into symmetric matrices
Cw2 = Cw2+triu(Cw2, 1)';
par.contactPar.weights = repelem(Cw1, ageBlockSizes, ageBlockSizes);        % expand blocks to create a 16x16 matrix that can multiplied elementwise with the contact matrix
par.contactPar.weightsChange = repelem(Cw2, ageBlockSizes, ageBlockSizes);

% Start date and time window for change of contact matrix - matrix will
% change linearly during the specified number of days followiong the start
% date
par.contactPar.changeDate = par.rampStart;      
par.contactPar.changeWindow = 70 + plusMinusInt(Theta.MRampDays, 20);
    



%---------------------------  VOC model ---------------------------------
par.vocWaneAmount = 0.39;                       % coefficient determining what fraction of each post-infection susceptibler compartment gets bumped down the immunity scale
par.vocWaneDate = datenum('20JUN2022');                       % date of VOC dominance
par.vocWaneWindow = 2;                                        % time window over which to apply VOC parameter changes



%------------------------- Immunity parameters --------------------------
par.waneRateMult = 1 + plusMinus(Theta.waneRate, 0.4);               % multiplier on waning rates
par.waneRate_RtoS = par.waneRateMult * 1/60*0.5;                      % rate of moving R to S
par.waneRate_StoS = par.waneRateMult * 0.009*0.5;                     % rate of moving from one S compartment to the next one with lower immunity

% Run Khoury/Golding submodel to generate immunity parameters
kLog = 2.94/log(10);               % steepness of logistic relationship betweem log titre and VE -  2.94/log(10) -  Khoury Table S5
no50_sympt = log(0.2);         %  determines mapping from titre to VE symptoms - log(0.2)  
no50_sev = log(0.03);            % determines mapping from titre to VE severe - log(0.04)  
no50_inf = log(0.3);            
no50_trans = log(1.1);          
logTitreRatio = log(10);               % ratio of titre from one compartment to next (don't change this without also changing waning rate and calibrating to Golding reulsts)

logTitre0_2 = log(0.2);                 % strength of immunity (measured as initial titre) after 2 doses
logTitre0_3 = log(0.4);                 % strength of immunity (measured as initial titre) after 3 doses
logTitre0_inf = log(0.8)+log(5);               % strength of immunity (measured as initial titre) after 0/1 doses + infection
logTitre0_inf_plus2 = log(3)+log(5);           % strength of immunity (measured as initial titre) after 2 doses + infection
logTitre0_inf_plus3 = log(7)+log(5);          % strength of immunity (measured as initial titre) after 3 doses + infection

minVEsev = 0.5;                 % immunity to hospitalisation and death (from vaccine or prior infection) cannot wane below this   

% set titre levels for each susceptible compartment 
logTitreSequence = logTitreRatio*[0, -1, -2, -3];
logTitreLevels = [-inf -inf    logTitre0_2 + logTitreSequence   logTitre0_3 + logTitreSequence     logTitre0_inf_plus3 + logTitreSequence]; 

% convert titre levels to immunity agaist each outcome
par.VEi = getKhouryVE(logTitreLevels, kLog, no50_sympt);
par.VEh = getKhouryVE(logTitreLevels, kLog, no50_sev);
par.VEt = zeros(1, parBase.nSusComp);
par.VEs = par.VEi;
par.VEf = par.VEh;

% apply minimum immunity constraint to severe outcomes
par.VEh(3:end) = max(minVEsev, par.VEh(3:end));
par.VEf(3:end) = max(minVEsev, par.VEf(3:end));

% Calculate the proporton of post-infection indiciduala with 0/1 or 2 doses
% who go to each of the 4 post-infection susceptible compartments
% this is done by solving an ODE to calculate the proportion of a
% post-infection cohort who are in each compartment at time t such that
% their average titre has dropped by the specified amount
logTitreDrop =  [logTitre0_inf, logTitre0_inf_plus2] - logTitre0_inf_plus3;
Y0 = getImmPars(logTitreRatio, logTitreDrop );
Y0 = Y0./sum(Y0, 2);

% Proportion going into W1, W2, W3, W4 post recovery
postRecovDist_Unvaxed = Y0(1, :);
postRecovDist_Vaxed2 = Y0(2, :);
postRecovDist_Vaxed3 = [1 0 0 0];

[par.waneNet_StoS, par.waneNet_RtoS, par.vaxNet] = getATCmatrices(postRecovDist_Unvaxed, postRecovDist_Vaxed2, postRecovDist_Vaxed3);


% ------------------------ VOC VE model ----------------------------------
% To model reduced VE for VOC set VOC_titreDrop < 1 
par.VOC_logTitreDrop = log(0.4);
logTitre0_2_VOC = par.VOC_logTitreDrop + logTitre0_2;                 
logTitre0_3_VOC = par.VOC_logTitreDrop + logTitre0_3;
logTitreLevels_VOC = [-inf -inf    logTitre0_2_VOC + logTitreSequence   logTitre0_3_VOC + logTitreSequence     logTitre0_inf_plus3 + logTitreSequence]; 
par.VEi_VOC = getKhouryVE(logTitreLevels_VOC, kLog, no50_sympt);
par.VEh_VOC = getKhouryVE(logTitreLevels_VOC, kLog, no50_sev);
par.VEt_VOC = zeros(1, parBase.nSusComp);
par.VEs_VOC = par.VEi_VOC;
par.VEf_VOC = par.VEh_VOC;

% apply minimum immunity constraint to severe outcomes
par.VEh_VOC(3:end) = max(minVEsev, par.VEh_VOC(3:end));
par.VEf_VOC(3:end) = max(minVEsev, par.VEf_VOC(3:end));




