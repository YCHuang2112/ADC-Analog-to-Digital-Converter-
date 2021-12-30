% clear;
addpath ('D:\LocalLaptop\õìëË\å§ãÜèä\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\õìëË\å§ãÜèä\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');

fsig = 13e6;  osr=16; n_sample=100*osr; fB=50e6*osr; fs=fB*2; ts=1/fs; 
% fvco=1e8; 
fvco=1e8*osr; 
tvco=1/fvco; vco_inverter_nums=33;
tsys_res=1/fvco;
vco_levels = (vco_inverter_nums-3)/2*2;
ideal_sqnr =  10*log10(3*(osr^3)*pi/(pi^3))+20*log10(vco_levels)+1.76;

t_sim = (n_sample-1)*ts;
%% Open Simulink diagram
% options=simset('RelTol', 1e-3, 'MaxStep', 1/fs/1000);
%% model selection for different Simulink versions
% open('vco.slx'); out = sim('vco.slx', t_sim, options); % 2015 Simulink
open('vco.slx'); out = sim('vco.slx', t_sim); % 2015 Simulink

t=0:ts:ts*(n_sample-1);
yt = out.vco_out';
yt = double(yt(end-n_sample+1:end));
% yt = yt - mean(yt);
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
