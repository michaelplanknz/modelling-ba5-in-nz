function [dates, doses1, doses2, doses3, doses4plus] = getVaccineData(fNamePart, vaxDataDate, par)

fName = "data/vaccine_data_" + fNamePart + ".csv";
fprintf('   Loading vaccination data:    %s\n', fName)
vaxData = readtable(fName);

fNameProj = "data/reshaped_b2_projections_final_2022-07-13.csv";
fprintf('    Loading 4th dose projections:     %s\n', fNameProj);
vaxProj = readtable(fNameProj);

% Join actuals and projections into one table
vaxMerged = outerjoin(vaxData, vaxProj, 'LeftKeys', 'dates', 'RightKeys', 'date', 'MergeKeys', true);

fprintf('4th dose actuals:\n')
fprintf(' %7d ', max(table2array(vaxMerged(:, 50:65))));
fprintf('\n')

[~, nCols] = size(vaxMerged);
for iCol = 2:nCols                  % overwrite NaNs by filling the last non-nan value down the column (skipping 1st column which is dates)
   ind = find(~isnan(table2array(vaxMerged(:, iCol))), 1, 'last');
   vaxMerged(ind+1:end, iCol) = vaxMerged(ind, iCol);
end



% For 4th doses overwrite actuals with projections where the latter is larger
vaxMerged.doses4_11 = max(vaxMerged.doses4_11, vaxMerged.cumul_doses_admin_50_54);
vaxMerged.doses4_12 = max(vaxMerged.doses4_12, vaxMerged.cumul_doses_admin_55_59);
vaxMerged.doses4_13 = max(vaxMerged.doses4_13, vaxMerged.cumul_doses_admin_60_64);
vaxMerged.doses4_14 = max(vaxMerged.doses4_14, vaxMerged.cumul_doses_admin_65_69);
vaxMerged.doses4_15 = max(vaxMerged.doses4_15, vaxMerged.cumul_doses_admin_70_74);
vaxMerged.doses4_16 = max(vaxMerged.doses4_16, vaxMerged.cumul_doses_admin_75_);

fprintf('4th dose projections:\n')
fprintf(' %7d ', max(table2array(vaxMerged(:, 50:65))));
fprintf('\n')



dates = datenum(vaxMerged.dates_date)';
doses1 =     [vaxMerged.doses1_1 vaxMerged.doses1_2  vaxMerged.doses1_3  vaxMerged.doses1_4  vaxMerged.doses1_5  vaxMerged.doses1_6  vaxMerged.doses1_7  vaxMerged.doses1_8  vaxMerged.doses1_9  vaxMerged.doses1_10  vaxMerged.doses1_11  vaxMerged.doses1_12  vaxMerged.doses1_13  vaxMerged.doses1_14  vaxMerged.doses1_15  vaxMerged.doses1_16 ];
doses2 =     [vaxMerged.doses2_1 vaxMerged.doses2_2  vaxMerged.doses2_3  vaxMerged.doses2_4  vaxMerged.doses2_5  vaxMerged.doses2_6  vaxMerged.doses2_7  vaxMerged.doses2_8  vaxMerged.doses2_9  vaxMerged.doses2_10  vaxMerged.doses2_11  vaxMerged.doses2_12  vaxMerged.doses2_13  vaxMerged.doses2_14  vaxMerged.doses2_15  vaxMerged.doses2_16 ];
doses3 =     [vaxMerged.doses3_1 vaxMerged.doses3_2  vaxMerged.doses3_3  vaxMerged.doses3_4  vaxMerged.doses3_5  vaxMerged.doses3_6  vaxMerged.doses3_7  vaxMerged.doses3_8  vaxMerged.doses3_9  vaxMerged.doses3_10  vaxMerged.doses3_11  vaxMerged.doses3_12  vaxMerged.doses3_13  vaxMerged.doses3_14  vaxMerged.doses3_15  vaxMerged.doses3_16 ];
doses4plus = [vaxMerged.doses4_1 vaxMerged.doses4_2  vaxMerged.doses4_3  vaxMerged.doses4_4  vaxMerged.doses4_5  vaxMerged.doses4_6  vaxMerged.doses4_7  vaxMerged.doses4_8  vaxMerged.doses4_9  vaxMerged.doses4_10  vaxMerged.doses4_11  vaxMerged.doses4_12  vaxMerged.doses4_13  vaxMerged.doses4_14  vaxMerged.doses4_15  vaxMerged.doses4_16 ];


dates = dates + par.vaccImmDelay;     % shift vaccination dates to allow for delay


iDate = find(datenum(dates) == par.date0);
nPad = par.tEnd - (length(dates) - iDate);
if nPad > 0
    dates = [dates, dates(end)+1:dates(end)+nPad];
    doses1 = [doses1(1:end, :); repmat(doses1(end, :), nPad, 1)];
    doses2 = [doses2(1:end, :); repmat(doses2(end, :), nPad, 1)];
    doses3 = [doses3(1:end, :); repmat(doses3(end, :), nPad, 1)];
    doses4plus = [doses4plus(1:end, :); repmat(doses4plus(end, :), nPad, 1)];
elseif nPad < 0
    dates = dates(1:end+nPad);
    doses1 = doses1(1:end+nPad, :);
    doses2 = doses2(1:end+nPad, :);
    doses3 = doses3(1:end+nPad, :);
    doses4plus = doses4plus(1:end+nPad, :);
end






dates = dates(iDate:end);
doses1 = doses1(iDate:end, :);
doses2 = doses2(iDate:end, :);
doses3 = doses3(iDate:end, :);
doses4plus = doses4plus(iDate:end, :);



