function [snrdB,ptotdB,psigdB,pnoisedB] = calcSNR(vout,f,fBL,fBH,w,N)
% SNR calculation in the time domain (P. Malcovati, S. Brigati)
%
% [snrdB,ptotdB] = calcSNR(vout,f,fBL,fBH,w,N)
% [snrdB,ptotdB,psigdB] = calcSNR(vout,f,fBL,fBH,w,N)
% [snrdB,ptotdB,psigdB,pnoisedB] = calcSNR(vout,f,fBL,fBH,w,N)
%
% vout:     	Sigma-Delta bitstream taken at the modulator output
% f:    		Normalized signal frequency (fs = 1)
% fBL:			Base-band lower limit frequency bins
% fBH:			Base-band upper limit frequency bins
% w:			Windowing vector
% N:			Number of samples
%
% snrdB: 	 	SNR in dB
% ptotdB:  		Sigma-Delta modulator output power spectral density (vector) in dB
% psigdB:	 	Extracted signal power spectral density (vector) in dB
% pnoisedB:		Noise power spectral density (vector) in dB

fBL=ceil(fBL);
fBH=ceil(fBH);
signal=(N/sum(w))*sinusx(vout(1:N).*w,f,N);	% Extracts sinusoidal signal
noise=vout(1:N)-signal;			            % Extracts noise components
%noise(1:2)=0; % added by Derek on 1/21/2011
stot=((abs(fft((vout(1:N).*w)'))).^2);		% Bitstream PSD
ssignal=(abs(fft((signal(1:N).*w)'))).^2;	% Signal PSD
snoise=(abs(fft((noise(1:N).*w)'))).^2;		% Noise PSD

pwsignal=sum(ssignal(fBL:fBH));	            % Signal power
pwnoise=sum(snoise(fBL:fBH));		        % Noise power

snr=pwsignal/pwnoise;
snrdB=dbp(snr);
norm=sum(stot(1:N/2))/sum(vout(1:N).^2)*N;	% PSD normalization
if nargout > 1
	ptot=stot/norm; 
	ptotdB=dbp(ptot);
end

if nargout > 2
	psig=ssignal/norm;
	psigdB=dbp(psig);
end

if nargout > 3
	pnoise=snoise/norm;
	pnoisedB=dbp(pnoise);
end
