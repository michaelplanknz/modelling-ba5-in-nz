function errorShadeFull(x, Y, qt, clr, ls, mkr)

% Plot median plus quantiles as shaded bands 

x = reshape(x, 1, numel(x));
[n, m] = size(Y);
nBands = length(qt)/2;
if nBands > 7
    fprintf(    'WARNING: errorShadeFull cannot plot >7 bands\n\n')
elseif  mod(nBands, 1) ~= 0
    fprintf(    'WARNING: errorShade expects qt to have an even number of elements')
end

clrCoeff = 0.75 - 0.2*(nBands/7)*linspace(-1, 1, nBands);

ci = quantile(Y, qt, 1);

for iBand = 1:nBands
    ind1 = iBand;
    ind2 = 2*nBands+1-iBand;
    inLowerFlag = ~isnan(ci(ind1, :));
    inUpperFlag = ~isnan(ci(ind2, :));

    xShade = [x(inLowerFlag), fliplr(x(inUpperFlag))];
    yShade = [ci(ind1, inLowerFlag), fliplr(ci(ind2, inUpperFlag))];
    fill(xShade, yShade, clrCoeff(iBand) + (1-clrCoeff(iBand))*clr, 'LineStyle', 'none', 'HandleVisibility', 'off', 'FaceAlpha', 0.5)
    hold on
end
plot(x, median(Y, 1), 'LineStyle', ls, 'Marker', mkr, 'Color', clr)
