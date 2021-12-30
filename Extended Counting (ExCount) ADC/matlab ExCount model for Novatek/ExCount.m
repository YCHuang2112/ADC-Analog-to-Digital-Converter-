%% Incremental ADC-  AC test
clear
addpath ('SDtoolbox');

%% Modulator Parameter
vref=0.9; vfs=vref*2; VFS=vfs; %vcm=vfs/2; vicm=vcm; % full-scale vfs=1.8V
vin_amp_dB=[-80:5:-10, -9.5:0.5:0];  
vin_amp=10.^(vin_amp_dB/20)*vref;

%% Frequency of Signal and Clock
fB=50e3; fnyq=2*fB; T_nyq=1/fnyq; MUP=1; %% nyquist output frequency
OSR=16; M=OSR;
fs=OSR*fnyq; Ts=1/fs;  fclk=fs; tclk=Ts;
fsig=13.0e3; Tsig=1/fsig; %% input signal frequency
n_sample=OSR*100; 
n_cycle=fsig * n_sample / fs;
n_nyq=n_cycle/fsig*fnyq;
t_sim= n_cycle/fsig; 


%% Open Simulink diagram
options=simset('RelTol', 1e-3, 'MaxStep', 1/fs);
%% model selection for different Simulink versions
% open('ExCount_naive_model_R2012b.mdl');
% sim('ExCount_naive_model_R2012b', t_sim, options); % 2012 Simulink
open('ExCount_naive_model.slx');
sim('ExCount_naive_model', t_sim, options); % 2019 Simulink

%% FFT @ Nyquist Rate
fbin_sig=fsig/fnyq; w=hann_pv(n_nyq); fbin_L=3; fbin_H=n_nyq*fB/fnyq;
fbin_sig_pre=fsig/fclk; w_pre=hann_pv(n_sample); fbin_L_pre=3; fbin_H_pre=n_sample*fB/fclk;
sqnrs=[];
ptots=[];
for m=1:length(vin_amp_dB)
    ADCout_n=ADCout(end-n_nyq+1:end,m);
    [sqnrs(end+1),ptots(:,end+1)]      = calcSNR(ADCout_n',fbin_sig,fbin_L,fbin_H,w,n_nyq);
end
    


%% peak SNDR v.s.OPAMP GAIN (Incremental ADC) 
peak_SQNR = max(sqnrs);
vin_dB_at_peak_SQNR_index = find(sqnrs==peak_SQNR);
vin_dB_at_peak_SQNR = vin_amp_dB(vin_dB_at_peak_SQNR_index(1));

figure(1); clf;
title (['SQNR v.s. Input Amplitude']);
plot(vin_amp_dB,sqnrs); 
text_handle= text(-70,90, sprintf('Peak SQNR = %4.1fdB\n@ Amplitude = %4.1fdB',peak_SQNR,vin_dB_at_peak_SQNR),'Color','b','FontSize',18);
xlabel ('Amplitude (dBFS)'); ylabel ('SQNR (dB)'); axis([-80 0 0 105]);

N_start=1;
N_stop=40;
figure(2); clf;
subplot(3,2,1);
plot(n1(N_start:N_stop,vin_dB_at_peak_SQNR_index));
legend('n1','FontSize',12);
subplot(3,2,2);
plot(n2(N_start:N_stop,vin_dB_at_peak_SQNR_index));
legend('n2','FontSize',12);
subplot(3,2,3);
plot(D1d(N_start:N_stop,vin_dB_at_peak_SQNR_index));
legend('D_1_d','FontSize',12);
subplot(3,2,4);
plot(D2(N_start:N_stop,vin_dB_at_peak_SQNR_index));
legend('D_2','FontSize',12);
subplot(3,2,5);
plot(Ds(N_start:N_stop,vin_dB_at_peak_SQNR_index));
legend('D_s','FontSize',12);
xlabel ('nth points','FontSize',14);
subplot(3,2,6);
plot(ADCout(N_start:N_stop,vin_dB_at_peak_SQNR_index));
legend('D','FontSize',12);
xlabel ('nth points','FontSize',14);
%% END