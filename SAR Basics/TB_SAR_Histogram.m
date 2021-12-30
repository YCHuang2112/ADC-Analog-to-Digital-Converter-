%% include Toolbox
% clear all;
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');
%% global parameter setting
% global voltage setting
vref=0.9;
vfs=2*vref;

%% User Parameters
% Cu_sigma_USER = 0.01; xlim_start=56; xlim_end=62; x_txt_loc=56.5;
% Cu_sigma_USER = 0.03; xlim_start=52; xlim_end=62; x_txt_loc=53;
Cu_sigma_USER = 0.05; xlim_start=42; xlim_end=62; x_txt_loc=44;
N_round_mismatch_USER = 100;



%% sar setting
% bit setting
% B > 1
B = 10;
if B <= 1
    error('B has to be larger than 1 (>1)#n');
end

%%% DAC setting
% Cp_array = [ 2.^[(B-2):-1:0 0] ];
% Cn_array = [ 2.^[(B-2):-1:0 0] ];
Cp_array = [140.253 70.2098 35.1033 17.55156 8.77391 4.3846 2.18671 1.10484 0.553204 0.550779];
Cn_array = [140.253 70.2098 35.1033 17.55156 8.77391 4.3846 2.18671 1.10484 0.553204 0.550779];
Cu_mu = 0.553204;       % Cu = Cu_mu (fF)
Cp = 0;                 % Parasitic Cap = Cp (fF)

Cu_sigma = Cu_sigma_USER; % Cap mismatch variance relative to Cu



C_nums = 2.^[(B-2):-1:0 0];
% for i = 1:length(C_nums)
%     Cp_array_mis(1,i) = sum(randn(1,C_nums(i)))*Cu_sigma + Cp_array(1,i);
%     Cn_array_mis(1,i) = sum(randn(1,C_nums(i)))*Cu_sigma + Cn_array(1,i);
% end
% 
% vdacp = [ Cp_array_mis(1:B-1) 0]/(sum(Cp_array_mis)+Cp)*vfs;
% vdacn = [ Cn_array_mis(1:B-1) 0]/(sum(Cn_array_mis)+Cp)*vfs;

%% input signal & sampling setting
fB=500e3; 
fnyq=2*fB; Tnyq=1/fnyq; %% nyquist output frequency
fclk = fnyq * B; Tclk = 1/fclk; 
N = 1000;

% input sine-signal setting
fsig = 173e3;
vin_amp_dB=-0.5; vin_amp = vfs*10^(vin_amp_dB/20); % single-ended amplitude
Cycle=fsig*N/fnyq;

if  floor(Cycle) ~= ceil(Cycle)
    error([sprintf("ERROR: (Cycle)=(%f) should be an integer\n",Cycle)]);
end
if gcd(N,Cycle) ~= 1
    error([sprintf("ERROR: (N,Cycle)=(%d,%d) should be relatively prime\n",N,Cycle)]);
end

%% simulation setting
% simulation setting
t_sim = Tnyq*N;
options=simset('RelTol', 1e-3, 'MaxStep', 1/fclk);



%%  Mismatch Histogram setting
N_round_mismatch = N_round_mismatch_USER;
snr_s = zeros(1,N_round_mismatch);

for i_round = 1:N_round_mismatch
    % DAC Mismatch setting
    for i = 1:length(C_nums)
        Cp_array_mis(1,i) = sum(randn(1,C_nums(i)))*Cu_sigma * Cu_mu + Cp_array(1,i);
        Cn_array_mis(1,i) = sum(randn(1,C_nums(i)))*Cu_sigma * Cu_mu + Cn_array(1,i);
    end

    vdacp = [ Cp_array_mis(1:B-1) 0]/(sum(Cp_array_mis)+Cp)*vfs;
    vdacn = [ Cn_array_mis(1:B-1) 0]/(sum(Cn_array_mis)+Cp)*vfs;

    % start simulation
    open('SAR.slx'); 
    out = sim('SAR', t_sim, options);

    % extract result
    sar_dout = out.SAR_DOUT(end-N+1:end,1)';

    % calculate snr
    fbin_sig=fsig/fnyq; w=hann_pv(N); fbin_L=3; fbin_H=N*fB/fnyq;
    [snr, ptot]      = calcSNR(sar_dout,fbin_sig,fbin_L,fbin_H,w,N);
    
    snr_s(i_round) = snr;
end

%% diagram representation
figure(1);
set(gcf, 'color', [1 1 1]);
hist(snr_s);
text(x_txt_loc,23, sprintf('mean  = %3.1fdB',mean(snr_s)),'Color','g','FontSize',36);
text(x_txt_loc,16, sprintf('sigma = %3.1fdB',sqrt(var(snr_s))),'Color','m','FontSize',36);
xlabel('SQNR [dB]','FontSize',20); ylabel('Count(s)','FontSize',20); 
% title({['Histogram of output SQNR'];['with \sigma = ',num2str(Cu_sigma*100),'% Mismatch (N=100)']}); 
title({['Histogram of output SQNR'];['under ',num2str(Cu_sigma*100),'%-\sigma',' Mismatch (N=100)']}); 
axis([xlim_start xlim_end 0 25]);
set(gca, 'FontSize',20);



% figure(1);
% plot(sar_dout);
% 
% figure(2); %clf; 
% set(gcf, 'color', [1 1 1]);
% hold on;
% plot(linspace(0,fnyq/2,N/2)/1e3, ptot(1:N/2), 'r','linewidth',2.0); 
% hold off;
% text_handle= text(10,-50, sprintf('SQNR = %4.1fdB',snr),'FontSize',32);
% xlabel('Frequency [kHz]','FontSize',20); ylabel('PSD [V^2/Hz]','FontSize',20);  axis([0 fnyq/2/1e3 -130 0]);
% title ({[ 'SAR 10-bit;  V_P_P= ',num2str(vin_amp_dB),'dBFS;  '],['F_s_i_g=',num2str(fsig/1e3),'kHz;  ','F_B=',num2str(fnyq/1e6/2),'MHz; ' ]});
% legend('SQNR');
% legend('boxoff');
% set(gca, 'FontSize',20);

