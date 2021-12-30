%% Feedforward Delta-Sigma Modulator only; 
%  not including decimation filter.  2nd-order & 3rd-order Mod.
clear;
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\MIT Design Tools\HSPICE_Toolbox\HspiceToolbox');

%% parameters
vref=0.9; vfs=vref*2; r1=8; OSR=16*r1; fB=50e3; fs=OSR*2*fB; Ts=1/fs;
M = OSR-1;
fnyq=2*fB; T_nyq=1/fnyq;
%n_sample=2^14;  n_cycle=63*8*2/1; 
%n_sample=2^14;  n_cycle=5*8*1/r1;
%fsig=(n_cycle/(n_sample/fs)); %vp_db=-10;  vp=vref*10^(vp_db/20);
fsig=21e3; n_sample=2^7*100; 
n_cycle=fsig * n_sample / fs;
n_nyq=n_cycle/fsig*fnyq;
%fsig=round(0.16/OSR*n_samle);
%bit_adc2=6; LSB2=vfs/(2^bit_adc2);
alpha2=pi/OSR/sqrt(3); %2nd order zero-optimization
g=2-2*cos(alpha2); g1=2-g; g2=1-g; % 2nd-order zero-optimization

%% Simulink Simulation
%% Circuit non-idealities parameters
%%%%%%%%%%%%%%%%%%% kT/C noise
k=1.38e-23;				% Boltzmann Constant
Temp=300;				% Absolute Temperature in Kelvin
Cs=0.01e-12;				% Integrating Capacitance of the first integrator;  1e-12=1p

%%%%%%%%%%%%%%%%%%% OpAmp Non-idealities
op_gain_db=60;   op_gain=10^(op_gain_db/20);     % Opamp DC gain
% alfa=(op_gain-1)/op_gain;		% A=Op-amp finite gain (alfa=(A-1)/A -> ideal op-amp alfa=1)
alfa = 1;
Amax = vref*1;          % Op-amp saturation value [V]
%sr=2000000e6;			% Op-amp slew rate [V/s]
%GBW=1500000e6;			% Op-amp GBW [Hz]
sr=2e6; GBW=1.5e6;
%%%%%%%%%%%%%%%%%%%
noise1=0;       		% 1st int. output noise std. dev. [V/sqrt(Hz)]
delta=0;        		% Random Sampling jitter (std. dev.) [s]
%%%%%%%%%%% misc settings
Ntransient=0; fclk=fs; Tclk=1/fclk;
index=20;  %% -6dB
vin_amp_dB=[-70:5:-10 -9:0.5:0];  vin_amp=10.^(vin_amp_dB(index)/20)*vref; %vin_amp is amplitude of input sine wave
t_sim= n_cycle/fsig + Ntransient*Ts;  %n_sim = n_cycle*T_sig/tclk + Ntransient;

%% Open Simulink diagram 
k49=1; kint2=1/1; k1=1/2.5; % adder & 9-lev Quantizer scaling;
mdac=[1 1 1 1 1 1 1 1]; % mdac is the coefficient to model the mismatch of the feedback DAC
%mdac=[1+0.099/100 1+0.075/100 1-0.068/100 1-0.11/100  1 1+0.085/100 1 1-0.095/100]*1.00;
options=simset('RelTol', 1e-3, 'MaxStep', 1/fclk);
open('CTIADC_mod2_bit1_ff.mdl'); sim('CTIADC_mod2_bit1_ff', t_sim, options);


%%
%%%% Output graphs
figure (3); % time-domain plot
subplot(4,1,1); plot(xin  ,'r'); title ('input');
subplot(4,1,2); plot(int1,'g'); title ('mod1 output');
subplot(4,1,3); plot(int2 ,'g'); title ('mod2 output');
subplot(4,1,4); plot(iadc_out,'b'); title ('output');


%%%%%%%%%%%%%%%%%%%%%%%%%%% FFT & Spectrum Graph %%%%%%%%%
%% FFT and graph with Malcovati Toolbox
fbin_sig=fsig/fnyq; w=hann_pv(n_nyq); fbin_L=3; fbin_H=n_nyq*fB/fnyq;
% f=fsig/fclk; bw=fclk/OSR/2;  w=hann_pv(n_sample); 
% for m=1:length(vin_amp_dB)    
    y2q=iadc_out(end-n_nyq+1:end);

    [snr_idc2q,ptot_idc2q]      = calcSNR(y2q',fbin_sig,fbin_L,fbin_H,w,n_nyq);
 
    
% end





%% Nyquist-rate PSD of Incremental ADC 
%index=1;
index=23; %-6dB
figure (21); %clf;  % frequency is at linear scale
plot(linspace(0,fnyq/2,n_nyq/2), ptot_idc2q(1:n_nyq/2), 'r','linewidth',1.0); hold on;
text_handle= text(floor(fB/8),-50, sprintf('SQNR = %4.1fdB',snr_idc2q),'FontSize',18);
xlabel('Frequency [Hz]'); ylabel('PSD [dB/Hz]'); title('Incremental ADC'); axis([0 fnyq/2 -130 0]);
title ([ 'IADC2;  V_P= ',num2str(vin_amp_dB(index)),'dBFS;  ','Freq=',num2str(fsig/1e3),'kHz;  ' ,'fs=',num2str(fs/1e6),'MHz; ' ]);
legend('SQNR');

% %% SNR v.s. Amplitude (Incremental ADC)
% figure (22); %clf; %% All IDC
% plot(vin_amp_dB,snr_idc2q,'--r*','linewidth',0.5); hold on;
% legend('SQNR', 'SNR');
% %title (['Peak idc2q SNR = ',num2str(snr_idc2q(index)),'  Overloaded Range = ',num2str(vin_amp_db(index)),'dB']);
% title ('SNR v.s. Amplitude');
% xlabel ('Amplitude (dBFS)'); ylabel ('SNR (dB)'); axis([-90 0 0 300]);
% legend('SQNR');


%% Histograms of the integrator outputs
nbins=length(int1);
[bin1,xx1]=histo(int1(:),    length(int1) );     
[bin2,xx2]=histo(int2(:),    length(int2) );  
% [bin3,xx3]=histo(int1_new(:,index),length(int1_new) ); 
%[bin3,xx3]=histo(int1_d,      length(int1_d) );  [bin4,xx4]=histo(int2_d,        length(int2_d) );
figure(8); clf;
subplot(2,1,1);plot(xx1, bin1);
grid on; title('1^s^t Integrator output'); xlabel('Voltage [V]'); ylabel('Occurrences'); xlim([-1.8 1.8]);
subplot(2,1,2), plot(xx2, bin2);
grid on; title('2^n^d Integrator output'); xlabel('Voltage [V]'); ylabel('Occurrences');  xlim([-1*vref vref]);
% subplot(2,2,3), plot(xx3, bin3);
% grid on; title('New AX49 1st Integrator output'); xlabel('Voltage step [V]'); ylabel('Occurrences'); xlim([-1*0.5 0.5]);
% grid on; title('New AX49 1st Integrator output'); xlabel('Voltage [V]'); ylabel('Occurrences'); xlim([-1.8 1.8]);
% subplot(2,2,4), plot(xx4, bin4);
% grid on; title('New AX49 2nd Integrator output'); xlabel('Voltage step [V]'); ylabel('Occurrences');  xlim([-1*vref vref]);
% grid on; title('New AX49 2nd Integrator output'); xlabel('Voltage [V]'); ylabel('Occurrences');  xlim([-1*vref vref]);


%% end

