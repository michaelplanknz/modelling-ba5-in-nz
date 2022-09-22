function Y = combineBands(X)

% Merges two bands into one, e.g. to convert 5 year age bands into 10 year
% age bands
% If X is an n x m x p array, Y = combineBands(X) will return a n x q x p 
% array where q = ceil(m/2) and Y(i,j,k) = X(i,2*j-1,k)+X(i,2*j,k) 

Y  = X(:, 1:2:end, :);
Y2 = X(:, 2:2:end, :);
nCols = size(Y2, 2);
Y(:, 1:nCols, :) = Y(:, 1:nCols, :) + Y2;

