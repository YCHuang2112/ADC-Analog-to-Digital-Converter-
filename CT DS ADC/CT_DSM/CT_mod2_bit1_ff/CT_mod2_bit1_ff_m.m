%% Feedforward Delta-Sigma Modulator only; 
%  not including decimation filter.  2nd-order & 3rd-order Mod.
clear;
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');
% addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\MIT Design Tools\HSPICE_Toolbox\HspiceToolbox');

%% parameters
vref=0.9; vfs=vref*2; r1=8; OSR=16*r1; fB=50e3; fs=OSR*2*fB; Ts=1/fs;
%n_sample=2^14;  n_cycle=63*8*2/1; 
%n_sample=2^14;  n_cycle=5*8*1/r1;
%fsig=(n_cycle/(n_sample/fs)); %vp_db=-10;  vp=vref*10^(vp_db/20);
fsig=21e3; n_sample=2^7*100; 
n_cycle=fsig * n_sample / fs;
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
alfa=(op_gain-1)/op_gain;		% A=Op-amp finite gain (alfa=(A-1)/A -> ideal op-amp alfa=1)
Amax = vref*1;          % Op-amp saturation value [V]
%sr=2000000e6;			% Op-amp slew rate [V/s]
%GBW=1500000e6;			% Op-amp GBW [Hz]
sr=2e6; GBW=1.5e6;
%%%%%%%%%%%%%%%%%%%
noise1=0;       		% 1st int. output noise std. dev. [V/sqrt(Hz)]
delta=0;        		% Random Sampling jitter (std. dev.) [s]
%%%%%%%%%%% misc settings
Ntransient=0; fclk=fs; tclk=1/fclk;
vin_amp_dB=[-70:5:-10 -9:0.5:0];  vin_amp=10.^(vin_amp_dB/20)*vref; %vin_amp is amplitude of input sine wave
t_sim= n_cycle/fsig + Ntransient*Ts;  %n_sim = n_cycle*T_sig/tclk + Ntransient;

%% Open Simulink diagram 
k49=1; kint2=1/1; k1=1/2.5; % adder & 9-lev Quantizer scaling;
mdac=[1 1 1 1 1 1 1 1]; % mdac is the coefficient to model the mismatch of the feedback DAC
%mdac=[1+0.099/100 1+0.075/100 1-0.068/100 1-0.11/100  1 1+0.085/100 1 1-0.095/100]*1.00;
options=simset('RelTol', 1e-3, 'MaxStep', 1/fclk);
open('CT_mod2_bit1_ff.mdl'); sim('CT_mod2_bit1_ff', t_sim, options);


%%
%%%% Output graphs
index=20;  %% -6dB
figure (3); % time-domain plot
subplot(4,1,1); plot(xin  (:,index),'r'); title ('input');
subplot(4,1,2); plot(int1 (:,index),'g'); title ('mod1 output');
subplot(4,1,3); plot(int2 (:,index),'g'); title ('mod2 output');
subplot(4,1,4); plot(y_ax49(:,index),'b'); title ('output');


%%%%%%%%%%%%%%%%%%%%%%%%%%% FFT & Spectrum Graph %%%%%%%%%
%% FFT and graph with Malcovati Toolbox
f=fsig/fclk; bw=fclk/OSR/2;  w=hann_pv(n_sample); 
for m=1:length(vin_amp_dB),
    %yy1=xin_adc(:,m); %low-distortion feedforward     
%     yy2=y_noise(:,m);  yy3=y_ax49(:,m); yy4=y_ax49_new(:,m);
    yy3=y_ax49(:,m);
%     fbin_L=1; fbin_H=n_sample*bw/fclk;  %%fbin_L=N_sample*1/fclk; fbin_H=N_sample*f_nyq/fclk;
    fbin_L=1; fbin_H=n_sample*bw/fclk;
%     ptot_noise(:,m)=zeros(1,n_sample);  ptot_ax49(:,m)=zeros(1,n_sample); ptot_ax49_new(:,m)=zeros(1,n_sample);
    ptot_ax49(:,m)=zeros(1,n_sample);
    
%     [snr_noise(m),    ptot_noise(:,m)]    = calcSNR(yy2',f,fbin_L,fbin_H,w,n_sample);
    [snr_ax49(m),     ptot_ax49(:,m)]     = calcSNR(yy3',f,fbin_L,fbin_H,w,n_sample);
%     [snr_ax49_new(m), ptot_ax49_new(:,m)] = calcSNR(yy4',f,fbin_L,fbin_H,w,n_sample);
    
end
% ptot_noise=ptot_noise+2.92; ptot_ax49=ptot_ax49+2.92;  ptot_ax49_new=ptot_ax49_new+2.92; %% normalization

%% Plot SNR figure;
figure(5);  clf; %% old AX49 v.s. new AX49
semilogx(linspace(0,fclk/2,n_sample/2),     ptot_ax49(1:n_sample/2,index), 'r'); hold on; %grid on;
% semilogx(linspace(0,fclk/2,n_sample/2), ptot_ax49_new(1:n_sample/2,index), 'b'); hold on;
%semilogx(linspace(0,fclk/2,n_sample/2), ptot_noise(1:n_sample/2,index), 'r'); 
title ( ['V_P= ',num2str(vin_amp_dB(index)),'dB;  ','Freq=',num2str(fsig/1e3),'kHz;  ','CLK=',num2str(fs/1e6),'MHz;  ','OSR=',num2str(OSR)] );
% legend('original AX49','revised AX49'); 
legend('original AX49');
xlabel('Frequency [Hz]'); ylabel('PSD [dB/Hz]');  axis([0 fclk/2 -180 0]); 
text_handle = text(floor(10e3),-30, sprintf('Original AX49 SQNR = %4.1fdB ',snr_ax49(index)));
% text_handle = text(floor(10e3),-50, sprintf('Revised AX49 SQNR = %4.1fdB ',snr_ax49_new(index)));

% figure(6);  clf;  % SQNR v.s. SNR with thermal noise
% semilogx(linspace(0,fclk/2,n_sample/2), ptot_ax49_new(1:n_sample/2,index), 'b'); hold on;
% semilogx(linspace(0,fclk/2,n_sample/2), ptot_noise(1:n_sample/2,index), 'r'); 
% title ( ['V_P= ',num2str(vin_amp_dB(index)),'dB;  ','Freq=',num2str(fsig/1e3),'kHz;  ','CLK=',num2str(fs/1e6),'MHz;  ','OSR=',num2str(OSR)] );
% legend('SQNR','SNR with thermal noise (inconsistent random seed)'); %legend('SQNR','SNR w/ thermal noise'); 
% xlabel('Frequency [Hz]'); ylabel('PSD [dB]');  axis([0 fclk/2 -140 0]); 
% text_handle = text(floor(10e3),-30, sprintf('AX49 SNR = %4.1fdB ',snr_noise(index)));
% text_handle = text(floor(10e3),-50, sprintf('AX49 SQNR = %4.1fdB ',snr_ax49_new(index)));
% s1=sprintf('   with thermal noise, SNR(dB)=%1.3f',snr_noise(index) );  disp(s1);


%% SNR v.s. Amplitude
figure(7);
plot(vin_amp_dB, snr_ax49,'bo-');  grid on; hold on;
% plot(vin_amp_dB, snr_ax49_new,'bo-');  grid on; hold on;
% plot(vin_amp_dB, snr_noise,'ro-'); 
% legend('SQNR','SNR');
title ('2nd-order modulator, SQNR (dB) v.s. Amplitude (dB)'); xlabel ('Amplitude (dB)'); ylabel ('SQNR (dB)');
ylim([0 100]);
s_peak=sprintf('@ %4.1fdBFS Amplitude = %4.1fdB\n',vin_amp_dB(index), snr_ax49(index)); text(-40,10, s_peak, 'fontsize',12);
% s_peak2=sprintf('SQNR @ %4.1fdBFS Amplitude = %4.1fdB\n',vin_amp_dB(index), snr_ax49_new(index)); text(-40,30, s_peak2, 'fontsize',12);

%% calculation of voltage steps
% for i=1:length(int1(:,index))-1
% int1_d(i) = int1_new(i+1,index)-int1_new(i,index);  
% end
% int1_d=int1_d';  
% % int2_d=int2_d';

%% Histograms of the integrator outputs
nbins=length(int1);
[bin1,xx1]=histo(int1(:,index),    length(int1) );     
[bin2,xx2]=histo(int2(:,index),    length(int2) );  
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

