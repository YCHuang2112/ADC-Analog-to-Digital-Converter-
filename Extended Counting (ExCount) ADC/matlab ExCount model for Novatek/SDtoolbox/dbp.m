function y = dbp(x)
% Calculates the input value in dB dbp(x) = 10*log10(x)
% (by S. Brigati, P. Malcovati)
%
% y = dbp(x)
%
% x:	Input
%
% y:	Output in dB

y = -Inf*ones(size(x));
nonzero = x~=0;
y(nonzero) = 10*log10(abs(x(nonzero)));
