function popCount = readPopnData(fs, nAgeGroups)

fprintf('   Loading population size data:    %s\n', fs)
popSizeData = readmatrix(fs); % Load NZ population structure from data folder

popCount = [popSizeData(1:nAgeGroups-1, 2); sum(popSizeData(nAgeGroups:end, 2))]; % Fill entries with population distribution % Aggregate 75+ age-groups
