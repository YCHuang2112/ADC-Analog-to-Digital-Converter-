clear all;
addpath ('D:\LocalLaptop\專題\研究所\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\專題\研究所\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');


filename = 'mod2_1bit_ff_vin_n6db_np_12800_fsig_21M.m';
% delimiterIn = '\t';
% A = importdata(filename, delimiterIn);
run(filename);
%  test = v__inp__
% test = test - 0.9
% % x = -pi:(2*pi/1710):pi;
% % test = sin(x)
%  test = dout;
 test = v__dout__(end-12800:end);

for i=1:length(test)
    if test(i)>0.9
        test(i)=1;
    else
        test(i)=-1;
    end
end

% vref=0.9; vfs=vref*2; VFS=vfs;
% vin_amp_db=-6;  
% vin_amp=10^(vin_amp_db/20)*vref;

fB=50e3; fnyq=2*fB; T_nyq=1/fnyq; MUP=1; %% nyquist output frequency
OSR=128; M=OSR-1;
fs=OSR*fnyq; Ts=1/fs;  
% fclk=fs; tclk=Ts;
fclk=fs*1; tclk=Ts;
fsig=21e3; Tsig=1/fsig; %% input signal frequency


% % % Ntransient=10;
% % % n_nyq=50; % n_nyq=n_cycle/fsig*fnyq;
%n_sample=OSR*n_nyq;
n_sample = length(test)-1;
% % % n_cycle=n_sample*fsig/fs;  
% % % t_sim= n_cycle/fsig + Ntransient*Ts;  n_sim=n_cycle*Tsig/Ts + Ntransient;
% % % mdb=0.1/100*0;  % double-sampling error
%%FB_index_in_freq = n_sample*fB/fclk;

fbin_sig_sp=fsig/fclk; w_sp=hann_pv(n_sample); fbin_L_sp=1; fbin_H_sp=n_sample*fB/fclk;
% m=1;
% yt=test(end-n_sample+1:end);
yt=test(2:end);
[snr_test, ptot_test] = calcSNR(yt, fbin_sig_sp, fbin_L_sp, fbin_H_sp, w_sp, n_sample);

figure(1); clf;
semilogx(linspace(0,fclk/2,n_sample/2), ptot_test(1:n_sample/2), 'b','linewidth',2.0);
% plot(linspace(0,fclk/2,n_sample/2), ptot_test(1:n_sample/2), 'b','linewidth',2.0);
text_handle= text(floor(fB/2),-60, sprintf('SQNR = %4.1fdB',snr_test),'Color','green');
xlabel('Frequency [Hz]'); ylabel('PSD [dB/Hz]'); xlim([0 fclk/2]); ylim([-140 0]);
title ([ 'mod2 1bit ADC;  V_P= ','-6','dBFS;  ','Freq=',num2str(fsig/1e3),'kHz;  ' ,'fs=',num2str(fs/1e6),'MHz; ','OSR= ', num2str(OSR) ]);

figure(2); clf;
semilogx(linspace(0,fclk/2,n_sample/2), ptot_test(1:n_sample/2), 'b','linewidth',2.0);
% plot(linspace(0,fclk/2,n_sample/2), ptot_test(1:n_sample/2), 'b','linewidth',2.0);
text_handle= text(floor(fB/2),-60, sprintf('SQNR = %4.1fdB',snr_test),'Color','green');
xlabel('Frequency [Hz]'); ylabel('PSD [dB/Hz]'); xlim([0 fclk/2/OSR]); ylim([-140 0]);
title ([ 'mod2 1bit ADC;  V_P= ','-6','dBFS;  ','Freq=',num2str(fsig/1e3),'kHz;  ' ,'fs=',num2str(fs/1e6),'MHz; ','OSR= ', num2str(OSR) ]);
%%

vout = yt;
f = fbin_sig_sp;
fBL = fbin_L_sp;
fBH = fbin_H_sp;
w = w_sp;
N = n_sample;

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