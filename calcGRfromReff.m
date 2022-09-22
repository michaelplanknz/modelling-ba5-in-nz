function r = calcGRfromReff(Reff, a, w)

% Calculate epidemic growth rate r from reproduction number Reff using
% Wallina-Lipsitch method with assumed generation time distribution w(a)
% where a is a vector of times (in days) and w is an associated vector of 
% probabilities  

w = w/trapz(a, w);
gMean = trapz(a, a.*w);
r = zeros(size(Reff));
myFn = @(x, Re)(Re - 1./trapz(a, w.*exp(-x.*a), 2 ));
xa = (-1:0.01:1)';
for ii = 1:length(Reff)
    r0 = (Reff(ii)-1)/gMean;    % use exponential GT solution for initialisation point
    r(ii) = fzero(@(x)myFn(x, Reff(ii)), r0);
end





