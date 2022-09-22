function [avgImmNPI, avgImmPI] = calcAvgImm(Stemp, R, VEi, par)

% Takes a n x m matrix Stemp whose n rows are time points and whose m = k x 16 columns correspond
% to the k categories (e.g. different susceptible compartments or groups of
% susceptibility classes) and 16 age groups.
% Columns are in order: [class1, age1], [class1, age2], ... [class1, age16], [class2, age1], ...  
% R is a n x 16 matrix for total recovered in each age class
%
% Returns two n x 16 matrices (one for non-previously infected and one for previously infected) containing the avg immunity in each age group at each time point.

[nDays, nCols] = size(Stemp);
nToSum = nCols/par.nAgeGroups;

AtimesVEi = Stemp.*repelem(VEi, 1, par.nAgeGroups);             % A * VEi
 
Ar = reshape(Stemp, nDays, par.nAgeGroups, nToSum);             % reshape to 3D arrays with dimensions time x age x susceptibility compartment
AtimesVEir = reshape(AtimesVEi, nDays, par.nAgeGroups, nToSum);


avgImmNPI = sum(AtimesVEir(:, :, 1:10), 3)./sum(Ar(:, :, 1:10), 3);     % summing over susceptible classes 1-10 (not previously infected)
avgImmPI =  min(1, max(0, (R + sum(AtimesVEir(:, :, 11:14), 3))./(R + sum(Ar(:, :, 11:14), 3) )));     % summing over susceptible classes 1-10 (not previously infected), denominator includes all in recovered compartments

