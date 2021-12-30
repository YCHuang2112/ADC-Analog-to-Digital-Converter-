%% Incremental ADC-  AC test
clear
addpath ('D:\LocalLaptop\õìëË\å§ãÜèä\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\õìëË\å§ãÜèä\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');

%% User Parameters
% opgain_dB_USER = 80;
opgain_dB_USER = [70:10:100];
SAR_GAIN_ERROR_ARR_USER = [0.94:0.01:1.06];

%% Modulator Parameter
%vfs=3.3; vref=vfs/2; vcm=vfs/2; vicm=0; %vicm=vcm;  % full-scale vfs=3.3V
vref=0.9; vfs=vref*2; VFS=vfs; %vcm=vfs/2; vicm=vcm; % full-scale vfs=2V
vin_amp_dB=[-80:5:-10, -9.5:0.5:0];  
vin_amp=10.^(vin_amp_dB/20)*vref;

%% Frequency of Signal and Clock
fB=50e3; fnyq=2*fB; T_nyq=1/fnyq; %% n
OSR=16; M=OSR;
fs=OSR*fnyq; Ts=1/fs;  fclk=fs; tclk=Ts;
fsig=21.0e3; Tsig=1/fsig; %% input signal frequency
n_sample=2^7*100; 
n_cycle=fsig * n_sample / fs;
n_nyq=n_cycle/fsig*fnyq;
t_sim= n_cycle/fsig;  n_sim=n_cycle*Tsig/Ts;

stage1_Qlev = 2;
orders = 2;
k = orders*2+1;
% ideal_snr = 20*log10(stage1_Qlev)+1.76+10*log10(k*OSR^k*pi/(pi^k));
%% voltage scaling;
a1=1/1;  a2=1/1;   % a1: input scaling factor; a2: scaling factor for adder
% a11=1/4; b11=4/1;  % IDC2 1-bit modulator coefficients
a11=1; b11=1;  % IDC2 1-bit modulator coefficients
a12=1;

%% error
k1=1/1; k2=1/1;  % k1&k2: adder scaling
mdac=[1 1 1 1  1 1 1 1];  % no DAC error
SAR_GAIN_ERROR_ARR = SAR_GAIN_ERROR_ARR_USER;

%% Circuit non-idealities parameters
adder_max=vref*1;     % saturation value at Adder
Amax=vref*1;          % Op-amp saturation value [V]
opgain_dB=opgain_dB_USER; opgain=10.^(opgain_dB/20); 


%% Open Simulink diagram
options=simset('RelTol', 1e-3, 'MaxStep', 1/fs);
%% model selection for different Simulink versions
open('ExCount_wi_SarGainError.mdl');

iadc_out_arr1 = []; 
for i=1:length(opgain)
    for j=1:length(SAR_GAIN_ERROR_ARR)
        A = opgain(i);
        SAR_GAIN_ERROR = SAR_GAIN_ERROR_ARR(j);
    %     alfa = alfa_arr(i);
        sim('ExCount_wi_SarGainError', t_sim, options); % 2015 Simulink
        iadc_out_arr1(i,j,:,:) = out_no_cal(end-n_nyq+1:end,:)';
    end
end

%% FFT @ Nyquist Rate
snr_iadc_out = zeros(length(opgain),length(SAR_GAIN_ERROR_ARR),length(vin_amp_dB));
ptot_iadc_out = zeros(length(opgain),length(SAR_GAIN_ERROR_ARR),length(vin_amp_dB),n_nyq);

fbin_sig=fsig/fnyq; w=hann_pv(n_nyq); fbin_L=3; fbin_H=n_nyq*fB/fnyq;
for i = 1:length(opgain)
    for j = 1:length(SAR_GAIN_ERROR_ARR)
        for k=1:length(vin_amp_dB)
            iadc_out1= reshape(iadc_out_arr1(i,j,k,:),[1,n_nyq]);

            [snr_iadc_out(i,j,k),ptot_iadc_out(i,j,k,:)]   = calcSNR(iadc_out1,fbin_sig,fbin_L,fbin_H,w,n_nyq);
        end
    end
end

%% SNR v.s. Amplitude (Incremental ADC)
line_type=["--r*","--y*","--g*","--b*","--m*"];

figure (1); clf; %%
set(gcf, 'color', [1 1 1]);
hold on;
% lengend_txt={};
for i_opgain = 1:length(opgain)
    out_size = size(snr_iadc_out);
    snr_iadc_out_for_specific_gain =  reshape(snr_iadc_out(i_opgain,:,:),[length(SAR_GAIN_ERROR_ARR),length(vin_amp)]);
    plot(SAR_GAIN_ERROR_ARR,max(snr_iadc_out_for_specific_gain,[],2)','linewidth',4,'DisplayName',sprintf("A_o = %4.1fdB",opgain_dB(i_opgain)) ); 
%     txt = sprintf("A_o=%4.1fdB",opgain_dB);
%     lengend_txt(end+1) = {txt};
end
% legend(lengend_txt);
legend('FontSize',14,'Location','northeast');
legend('boxoff');
title ('SQNR  V.S.  SAR\_GAIN\_ERROR \alpha','FontSize',32);
xlabel ('SAR\_GAIN\_ERROR \alpha','FontSize',32); ylabel ('SQNR (dB)','FontSize',32); 
set(gca, 'FontSize',20);
hold off;

%% END