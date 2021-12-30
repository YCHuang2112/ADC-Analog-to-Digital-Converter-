clear;
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');

fsig = 13e6; n_sample=100; fB=50e6; fs=fB*2; ts=1/fs; 
Ts = ts;
vref=1.8; vfs=vref*2;
Vsig_dB = 0;
vin_amp=vref*10^(-Vsig_dB/20);


t_sim = (n_sample-1+20)*ts;
%% Open Simulink diagram
% options=simset('RelTol', 1e-3, 'MaxStep', 1/fs/1000);
%% model selection for different Simulink versions
% open('vco.slx'); out = sim('vco.slx', t_sim, options); % 2015 Simulink
open('pipeline.mdl'); out = sim('pipeline.mdl', t_sim); % 2015 Simulink

t=0:ts:ts*(n_sample-1);
yt = out.pipeline_out(end-n_sample+1:end)';
yt = yt - mean(yt);
fbin_sig_sp=fsig/fs; w_sp=hann_pv(n_sample); fbin_L_sp=3; fbin_H_sp=n_sample*fB/fs;
[snr_test, ptot_test] = calcSNR(yt, fbin_sig_sp, fbin_L_sp, fbin_H_sp, w_sp, n_sample);
enob = (snr_test-1.76)/6.02;

figure(1);
set(gcf, 'color', [1 1 1]);
plot(0:ts:(n_sample-1)*ts,yt);
figure(2);
set(gcf, 'color', [1 1 1]);
plot( linspace(0,fs/2/1e6,n_sample/2+1),ptot_test(1:n_sample/2+1));
text_handle= text(floor(fB/2.5/1e6),-20, sprintf('SQNR = %4.1f dB',snr_test),'Color','green','FontSize',24);
xlabel('Frequency [MHz]','FontSize',24); ylabel('PSD [dB/Hz]','FontSize',24);  axis([0 fB/1e6 -100 0]); 
title ({ ['3-Qlev 11-Stage Pipeline ADC',''];['V_P_P = ',num2str(Vsig_dB),' dBFS;  ','F_s_i_g = ',num2str(fsig/1e6),' MHz;  ' ,'F_B = ',num2str(fs/2/1e6),' MHz; ']},'FontSize',14);
set(gca, 'FontSize',13);
set(gcf, 'Position',  [488   342   580   420]);
