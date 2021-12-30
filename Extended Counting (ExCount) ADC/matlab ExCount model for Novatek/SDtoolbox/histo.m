function [bin,xx] = histo(y,N)
% Bins the elements of the input vector Y into N equally spaced containers
% within the range of Y (by S. Brigrati, P. Malcovati)
%
% [bin,xx] = histo(y,N)
%
% y:			Input vector
% N:			Number of bins
%
% bin:			Vector containing the number of occurrences for each bin
% xx:			Vector containing the position of each bin

range=ceil(2*max(y)*10)/10;
dx = range/N;
for i=1:N
	x(i) = -range/2 + (i-1)*dx + dx/2;
end
[bin,xx]=hist(y,x);
