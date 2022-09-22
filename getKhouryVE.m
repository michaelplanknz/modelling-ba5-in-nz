function VE = getKhouryVE(logTitreLevels, k, no50)

% Calculate immunity from log titre according to the Khoury et al model

VE = 1./(1+ exp(-k*(logTitreLevels-no50)));

