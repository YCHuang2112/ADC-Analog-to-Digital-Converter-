%% Incremental ADC-  AC test
clear
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');

%% Modulator Parameter
%vfs=3.3; vref=vfs/2; vcm=vfs/2; vicm=0; %vicm=vcm;  % full-scale vfs=3.3V
vref=0.9; vfs=vref*2; VFS=vfs; %vcm=vfs/2; vicm=vcm; % full-scale vfs=2V
dc_select=0; vin_dc=0.4; % 0=AC input; 1=DC input
%vin_amp_db=-6;  % for sine ac input only.
vin_amp_db=[-80:5:-10, -9.5:0.5:0];  
index=23; %-6dB
vin_amp=10.^(vin_amp_db(index)/20)*vref;

%% Frequency of Signal and Clock
fB=50e3; fnyq=2*fB; T_nyq=1/fnyq; MUP=1; %% nyquist output frequency
%OSR=400/4; M=OSR-1; %M1=64; M2=OSR-1-M1;
OSR=64*2; M=OSR-1;
fs=OSR*fnyq; Ts=1/fs;  fclk=fs; Tclk=Ts;
fsig=21.0e3; Tsig=1/fsig; Ntransient=0;%% input signal frequency
n_sample=2^7*100; 
n_cycle=fsig * n_sample / fs;
% n_cycle=fsig/1e3*4;  n_nyq=n_cycle/fsig*fnyq;
n_nyq=n_cycle/fsig*fnyq;
t_sim= n_cycle/fsig + Ntransient*Ts;  n_sim=n_cycle*Tsig/Ts + Ntransient;
mdb=0.1/100*0;  % double-sampling error

stage1_Qlev = 2;
orders = 4;
k = orders*2+1;
ideal_snr = 20*log10(stage1_Qlev)+1.76+10*log10(k*OSR^k*pi/(pi^k));
%% voltage scaling;
a1=1/1;  a2=1/1;   % a1: input scaling factor; a2: scaling factor for adder
% a11=1/3; b11=1/1;  % IDC2 1-bit modulator coefficients
a11=1/4; b11=4/1;  % IDC2 1-bit modulator coefficients
aff=1/1; % bff=[1/2 2/2 1/8];  %IDC1 feedforward scaling

%% error
% k1=1/3; k2=1/1;  % k1&k2: adder scaling a5 9-lev and 5-lev respectively; 
k1=1/1; k2=1/1;  % k1&k2: adder scaling a5 9-lev and 5-lev respectively; 
mdac=[1 1 1 1  1 1 1 1];  % no DAC error
mdac3=[1 1]; %mdac3=1+[0 1.051]/1000;
err_dac=1.051/1000*0; % tri-level DAC error

%% Circuit non-idealities parameters
adder_max=vref*1;     % saturation value at Adder
Amax=vref*1;          % Op-amp saturation value [V]
k=1.38e-23;				% Boltzmann Constant
Temp=300;				% Absolute Temperature in Kelvin
Cs=33.55e-15;			% Integrating Capacitance of the first integrator
opgain_dB=86; opgain=10^(opgain_dB/20);
alfa=(opgain-1)/opgain;	% opgain=Op-amp finite gain (alfa=(A-1)/A -> ideal op-amp alfa=1)
% alfa=1;         	% opgain=Op-amp finite gain (alfa=(A-1)/A -> ideal op-amp alfa=1)
sr=2000000e6;			% Op-amp slew rate [V/s]
GBW=1.5e6;			    % Op-amp GBW [Hz]
noise1=0;       		% 1st int. output noise std. dev. [V/sqrt(Hz)]
delta=1e-15;        	% Random Sampling jitter (std. dev.) [s]

%% Open Simulink diagram
options=simset('RelTol', 1e-3, 'MaxStep', 1/fs);
%% model selection for different Simulink versions
open('Two_Stage_DTIADC_CTIADC.mdl'); sim('Two_Stage_DTIADC_CTIADC', t_sim, options); % 2015 Simulink

%% FFT @ Nyquist Rate
fbin_sig=fsig/fnyq; w=hann_pv(n_nyq); fbin_L=3; fbin_H=n_nyq*fB/fnyq;
fbin_sig_pre=fsig/fclk; w_pre=hann_pv(n_sample); fbin_L_pre=3; fbin_H_pre=n_sample*fB/fclk;
% for m=1:length(vin_amp_db)
y2q=idc2q(end-n_nyq+1:end);
y2q_pre=idc2q_pre(end-n_sample+1:end);

[snr_idc2q,ptot_idc2q]      = calcSNR(y2q',fbin_sig,fbin_L,fbin_H,w,n_nyq);

[snr_idc2q_pre,ptot_idc2q_pre]      = calcSNR(y2q_pre',fbin_sig_pre,fbin_L_pre,fbin_H_pre,w_pre,n_sample);
% end
%ptot_idc2q=ptot_idc2q+5.28; ptot_idc2=ptot_idc2+5.28; %ptot_idc2dwa=ptot_idc2dwa+5.28; %ptot_idc2dsm=ptot_idc2dsm+4.48; 


% %% PSD of DS Modulator
% %index=1;
% index=23; %-6dB
% figure (11); %clf;  % frequency is at linear scale
% % plot(linspace(0,fclk/2,n_sample/2), ptot_idc2q_pre(1:n_sample/2,index), 'r','linewidth',1.0); hold on;
% % plot(linspace(0,fclk/2,n_sample/2), ptot_idc2_pre(1:n_sample/2,index), 'b','linewidth',2.0); hold on; axis([0 fclk/2 -130 0]);
% plot(linspace(0,fnyq/2,n_nyq/2), ptot_idc2q_pre(1:n_nyq/2,index), 'r','linewidth',1.0); hold on;
% plot(fB, 'b','linewidth',1.0)
% axis([0 fnyq/2 -130 0]);
% text_handle= text(floor(fB/8),-40, sprintf('SQNR = %4.1fdB',snr_idc2q_pre(index)),'FontSize',18);
% xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]'); title('Incremental ADC preprocess'); 
% title ([ 'IADC2;  V_P= ',num2str(vin_amp_db(index)),'dBFS;  ','Freq=',num2str(fsig/1e3),'kHz;  ' ,'fs=',num2str(fs/1e6),'MHz; ' ]);
% legend('SQNR pre');
% 
% %% SNR v.s. Amplitude (DS Modulator)
% figure (12); %clf; %% All IDC
% plot(vin_amp_db,snr_idc2q_pre,'--r*','linewidth',0.5); hold on;
% legend('SQNR pre');
% %title (['Peak idc2q SNR = ',num2str(snr_idc2q(index)),'  Overloaded Range = ',num2str(vin_amp_db(index)),'dB']);
% title ('SNR v.s. Amplitude');
% xlabel ('Amplitude (dBFS)'); ylabel ('SNR (dB)'); axis([-90 0 0 90]);
% legend('SQNR pre');




%% Nyquist-rate PSD of Incremental ADC 
%index=1;
% index=23; %-6dB
figure (21); %clf;  % frequency is at linear scale
plot(linspace(0,fnyq/2,n_nyq/2), ptot_idc2q(1:n_nyq/2), 'r','linewidth',1.0); hold on;
text_handle= text(floor(fB/8),-50, sprintf('SQNR = %4.1fdB',snr_idc2q),'FontSize',18);
xlabel('Frequency [Hz]'); ylabel('PSD [V^2/Hz]'); title('Incremental ADC'); axis([0 fnyq/2 -160 0]);
title ([ 'IADC2;  V_P= ',num2str(vin_amp_db(index)),'dBFS;  ','Freq=',num2str(fsig/1e3),'kHz;  ' ,'fs=',num2str(fs/1e6),'MHz; ' ]);
legend('SQNR');

% %% SNR v.s. Amplitude (Incremental ADC)
% figure (22); %clf; %% All IDC
% plot(vin_amp_db,snr_idc2q,'--r*','linewidth',0.5); hold on;
% legend('SQNR', 'SNR');
% %title (['Peak idc2q SNR = ',num2str(snr_idc2q(index)),'  Overloaded Range = ',num2str(vin_amp_db(index)),'dB']);
% title ('SNR v.s. Amplitude');
% xlabel ('Amplitude (dBFS)'); ylabel ('SNR (dB)'); axis([-90 0 0 150]);
% legend('SQNR');


%% DS Modulator histogram
figure(3); clf; % voltage swing + step
hist_x=[-1*vref:0.005:vref];
subplot(6,1,1), hist(int1q(:),hist_x); %plot(xx1, bin1);
grid on; title('INT1 histogram'); xlabel('Voltage [V]'); ylabel('Occurrences');xlim([-1*vref vref]);
subplot(6,1,2), hist(int2q(:),hist_x);
grid on; title('INT2 histogram'); xlabel('Voltage [V]'); ylabel('Occurrences');xlim([-1*vref vref]);
subplot(6,1,3), hist(adderq(:),hist_x);
grid on; title('ADDER histogram'); xlabel('Voltage [V]'); ylabel('Occurrences');xlim([-1*vref vref]);

subplot(6,1,4), hist(s2_int1q(:),hist_x); %plot(xx1, bin1);
grid on; title('S2\_INT1 histogram'); xlabel('Voltage [V]'); ylabel('Occurrences');xlim([-1*vref vref]);
subplot(6,1,5), hist(s2_int2q(:),hist_x);
grid on; title('S2\_INT2 histogram'); xlabel('Voltage [V]'); ylabel('Occurrences');xlim([-1*vref vref]);
subplot(6,1,6), hist(s2_adderq(:),hist_x);
grid on; title('S2\_ADDER histogram'); xlabel('Voltage [V]'); ylabel('Occurrences');xlim([-1*vref vref]);


%% text display
% s2q =sprintf('idc2 peak SQNR=%1.2f',snr_idc2q(index));  disp(s2q);


%% END