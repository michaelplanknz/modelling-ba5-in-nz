function [C_detBal, popDistBench] = getC_detBal(nAgeGroups)

% Construct the modified version of the Prem et al matrix that satisfies
% the detailed balance equations for the specified age distribution, as
% described in Steyn et al (2021)

%------------- Load Contact Matrix and Define NGM --------------
fs = 'data/nzcontmatrix.xlsx';
fprintf('   Loading contact matrix:    %s\n', fs)
C = readmatrix(fs); % Get Prem et al contact matrix from data folder



fs = 'data/popsize_national.xlsx';      % This should *ALWAYS* be the national population distribution 'popsize_national.xlsx'
fprintf('Reading benchmark (national) population size data:\n') 
popCountBench = readPopnData(fs, nAgeGroups);
popDistBench = popCountBench/sum(popCountBench);

C_detBal = 0.5*(C + (popDistBench.')./(popDistBench) .* (C.') );




