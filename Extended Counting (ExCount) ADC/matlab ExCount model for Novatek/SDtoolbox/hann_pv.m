function w = hann_pv(n)
% Hanning window for FFT (by S. Brigati, P. Malcovati)
%
% w = hann_pv(n)
%
% x:	Length of the window (number of Points)
%
% w:	Hanning window of length n

w = .5*(1 - cos(2*pi*(0:n-1)/n) );
