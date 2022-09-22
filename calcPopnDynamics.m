function dYdt = calcPopnDynamics(Y, b, Mu, r)

% calculate the contribution of population dynamics (births b, deaths Mu and
% ageing r) to the rate of change in each age group 

[~, nCols] = size(Y);

dYdt = [b; r*Y(1:end-1, :)] - [r*Y(1:end-1, :); zeros(1, nCols)] - Mu.*Y;
