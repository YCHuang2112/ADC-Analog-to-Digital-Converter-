clear;
addpath ('D:\LocalLaptop\õìëË\å§ãÜèä\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\õìëË\å§ãÜèä\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');

fsig = 1.3e6; n_sample=100; fB=5e6; osr =1; fs=fB*2; ts=1/fs; 
fvco = 1000e6; fsens=1000e6; 
tsys_res=1e-11;

t_sim = (n_sample-1)*ts;
%% Open Simulink diagram
% options=simset('RelTol', 1e-3, 'MaxStep', 1/fs/1000);
%% model selection for different Simulink versions
% open('vco.slx'); out = sim('vco.slx', t_sim, options); % 2015 Simulink
open('vco_counter.slx'); out = sim('vco_counter.slx', t_sim); % 2015 Simulink

t=0:ts:ts*(n_sample-1);
yt = out.vco_out';
yt = yt - mean(yt);
fbin_sig_sp=fsig/fs; w_sp=hann_pv(n_sample); fbin_L_sp=3; fbin_H_sp=n_sample*fB/fs/osr;
[snr_test, ptot_test] = calcSNR(yt, fbin_sig_sp, fbin_L_sp, fbin_H_sp, w_sp, n_sample);
enob = (snr_test-1.76)/6.02;

figure(1);
set(gcf, 'color', [1 1 1]);
plot(0:ts:t_sim,yt);
figure(2);
set(gcf, 'color', [1 1 1]);
plot( 0:fs/n_sample:fs-fs/n_sample,ptot_test);
text_handle= text(floor(fB/20),20, sprintf('SQNR = %4.1fdB',snr_test),'Color','green','FontSize',24);
