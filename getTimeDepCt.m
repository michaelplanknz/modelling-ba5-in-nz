function Ct = getTimeDepCt(t, dailyCaseRate, par)


% Time-dependent control function
Cti = interp1( par.date0 + [0:1:par.tEnd], par.Ct, t);

% Dynamic response function (0 = under trigger point (or before start date), 1 = over trigger)
dynResp = normcdf( dailyCaseRate - par.responseCaseTrig, 0, par.responseTrigSmoothWindow*par.responseCaseTrig) .* normcdf( t-par.responseStartDate, 0, par.responseStartWindow);
Ct = (1-dynResp).*Cti + dynResp*par.responseCt;

