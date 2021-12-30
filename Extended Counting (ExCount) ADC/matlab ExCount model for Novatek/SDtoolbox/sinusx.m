function outx = sinusx(in,f,n)
% Extracts of a sinusoidal signal (S. Brigati, P. Malcovati)
%
% outx = sinusx(in,f,n)
%
% in:		Input data vector
% f:		Normalized input signal frequency
% n:		Number of simulation points
%
% outx:		Sinusoidal signal

sinx=sin(2*pi*f*[1:n]);
cosx=cos(2*pi*f*[1:n]);
in=in(1:n);
a1=2*sinx.*in;
a=sum(a1)/n;
b1=2*cosx.*in;
b=sum(b1)/n;
outx=a.*sinx + b.*cosx;