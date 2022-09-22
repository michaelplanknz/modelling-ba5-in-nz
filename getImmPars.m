function Y0 = getImmPars(logTitreRatio, logTitreDrop)

% Calculate probabilities of entering each of the four post-infection
% susceptible parameters for people who are either 0/1-dose vaccinated or
% 2-dose vaccinated, as specified by their log titre relative to someone
% who is 3-dose vaccinated

w = 0.009;           % Output should be insensitve to this parameter as it is only used to find the solution Y* of the transition ODE at a point in time which the average titre is a given level
                     % The time taken to reach this level will depend directly on w but the solutiopn Y* at that time shouldn't

ic = [1; 0; 0; 0];
tSpan = [0 360];
myRHS = @(t, y)( [-w*y(1); w*(y(1)-y(2)); w*(y(2)-y(3)); w*y(3)] );     % solve cohort ODE
[t, Y] = ode45(myRHS, tSpan, ic);

logTitreSequence = logTitreRatio*[0, -1, -2, -3];

avgRelLogTitre = sum(Y.*logTitreSequence, 2);        % average titre at time t relative to time 0

% Find time at which titre for an infected+boosted person decays to
% expected initial titre for an infected+(doulbe vaxed or unvaxed) person,
% as specified by the relative values in the input variable titreDrop
% tOffset is a vector of such times for each value specified in titreDrop
tOffset       = interp1(avgRelLogTitre, t, logTitreDrop);

% State of ODE at those times - can use later to reinitialise ODE for these
% classes of immunity
% Each row of Y0 is a vector of 4 values for each value specified in
% titreDrop
Y0 = interp1(t, Y, tOffset);





