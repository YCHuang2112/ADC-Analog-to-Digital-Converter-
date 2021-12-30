clear;
%% include Toolbox
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\›“‘è\Œ¤‹†Š\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');
%% global parameter setting
% global voltage setting
vref=0.9;
vfs=2*vref;

%% sar setting
% bit setting
% B > 1
B = 8; %% B == # Cycles for one Output Conversion
if B <= 1
    error('B has to be larger than 1 (>1)#n');
end

%%% DAC setting
Cp_array = [ 2.^[(B-2):-1:0 0] ];
Cn_array = [ 2.^[(B-2):-1:0 0] ];
% Cp_array = [140.253 70.2098 35.1033 17.55156 8.77391 4.3846 2.18671 1.10484 0.553204 0.550779];
% Cn_array = [140.253 70.2098 35.1033 17.55156 8.77391 4.3846 2.18671 1.10484 0.553204 0.550779];
Cu_mu = 0.553204;       % Cu = Cu_mu (fF)
Cp = 0;                 % Parasitic Cap = Cp (fF)

Cu_sigma = 0.00 ; % Cap mismatch variance relative to Cu
C_nums = 2.^[(B-2):-1:0 0];
for i = 1:length(C_nums)
    Cp_array(1,i) = sum(randn(1,C_nums(i)))*Cu_sigma * Cu_mu + Cp_array(1,i);
    Cn_array(1,i) = sum(randn(1,C_nums(i)))*Cu_sigma * Cu_mu + Cn_array(1,i);
end

vdacp = [ Cp_array(1:B-1) 0]/(sum(Cp_array)+Cp)*vref;
vdacn = [ Cn_array(1:B-1) 0]/(sum(Cn_array)+Cp)*vref;

%% input signal & sampling setting
fB=50e3; OSR=8;
fnyq=2*fB*OSR; Tnyq=1/fnyq; %% nyquist output frequency
fclk=fnyq*B; Tclk=1/fclk; 
N=100*OSR; 

% input sine-signal setting
fsig = 23e3;
vin_amp_dB=-0.5; vin_amp = vref*10^(vin_amp_dB/20); % single-ended amplitude
Cycle=fsig*N/fnyq;

if  fsig > fB
    error([sprintf("ERROR: Fsig (%f) > FB (%f). Fsig should be within the BandWidth\n",fsig,fB)]);
end
if  floor(Cycle) ~= ceil(Cycle)
    error([sprintf("ERROR: (Cycle)=(%f) should be an integer\n",Cycle)]);
end
if gcd(N,Cycle) ~= 1
    error([sprintf("ERROR: (N,Cycle)=(%d,%d) should be relatively prime\n",N,Cycle)]);
end

%% simulation setting
% simulation setting
t_sim = Tnyq*N;
options=simset('RelTol', 1e-3, 'MaxStep', 1/fclk/10);

% start simulation
open('SAR_NS_v1.slx'); 
out = sim('SAR_NS_v1', t_sim, options);

%% result analyzing
% extract result
sar_dout = out.SAR_DOUT(end-N+1:end,1)';

% calculate snr
fbin_sig=fsig/fnyq; w=hann_pv(N); fbin_L=3; fbin_H=N*fB/fnyq;
[snr, ptot]      = calcSNR(sar_dout,fbin_sig,fbin_L,fbin_H,w,N);

%% diagram representation
figure(1);
plot(sar_dout);

figure(2); %clf; 
hold on;
% plot(linspace(0,fnyq/2,N/2+1), ptot(1:N/2+1), 'r','linewidth',1.0); 
semilogx(linspace(0,fnyq/2,N/2+1), ptot(1:N/2+1), 'r','linewidth',1.0); 
ptot(1:fbin_L-1)=-200;
Number_sig_bins = 2;
loc_sig_bin = fbin_sig * N + 1;
% plot(linspace(0,fnyq/2,N/2+1), 10*log10(cumsum([10.^(ptot(1:loc_sig_bin-Number_sig_bins)/10); zeros((Number_sig_bins-1)*2+1,1); 10.^(ptot(loc_sig_bin+Number_sig_bins:N/2+1)/10); ])), 'g','linewidth',2.0);
semilogx(linspace(0,fnyq/2,N/2+1), 10*log10(cumsum([10.^(ptot(1:loc_sig_bin-Number_sig_bins)/10); zeros((Number_sig_bins-1)*2+1,1); 10.^(ptot(loc_sig_bin+Number_sig_bins:N/2+1)/10); ])), 'g','linewidth',2.0);
text_handle= text(floor(fB/2),-50, sprintf('SQNR = %4.1fdB',snr),'FontSize',18);
xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]'); title('Incremental ADC'); axis([0 fnyq/2 -130 0]);
title ([ 'SAR;  V_P= ',num2str(vin_amp_dB),'dBFS;  ','Fsig=',num2str(fsig/1e3),'kHz;  ' ,'Conversion Rate (Fconv)=',num2str(fnyq/1e6),'MHz; ' ]);
legend('SQNR');
hold off;