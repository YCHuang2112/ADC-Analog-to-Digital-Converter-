%% include Toolbox
addpath ('D:\LocalLaptop\����\������\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\����\������\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');

%% sar setting
% global voltage setting
vfs=1.8;

% bit setting
% B > 1
B = 4;
if B <= 1
    error('B has to be larger than 1 (>1)#n');
end

%%% DAC setting
% 1. canonical SAR  cdac's equivalent voltage dac
% vdacp = [2.^(-1:-1:-B)]*vfs; -- for DAC P
% vdacn = [2.^(-1:-1:-B)]*vfs; -- for DAC N

% 2. canonical SAR  cdac's to voltage dac (vdac) trasnfer
Cp_array = 2.^[(B-1):-1:0 0];
Cn_array = 2.^[(B-1):-1:0 0];
vdacp = Cp_array(1:B)/sum(Cp_array)*vfs;
vdacn = Cn_array(1:B)/sum(Cn_array)*vfs;

% 3. SAR  cdac's to vdac with mismatch
% Cu_mu = 1;       % Cu = Cu_mu (fF)
% Cu_sigma = 0.01; % Cap mismatch variance relative to Cu
% 
% C_nums = 2.^[(B-1):-1:0 0];
% Cp_array = zeros(size(C_nums));
% Cn_array = Cp_array;
% 
% for i = 1:length(C_nums)
%     Cp_array(1,i) = sum(randn(1,C_nums(i)))*Cu_sigma + Cu_mu*C_nums(i);
%     Cn_array(1,i) = sum(randn(1,C_nums(i)))*Cu_sigma + Cu_mu*C_nums(i);
% end
% vdacp = Cp_array(1:B)/sum(Cp_array)*vfs;
% vdacn = Cn_array(1:B)/sum(Cn_array)*vfs;

% 4. MCS (Top Plate VCM-based with using half-Clsb)
% Cp_array = [ 2.^[(B-3):-1:0 0] ];
% Cn_array = [ 2.^[(B-3):-1:0 0] ];
% vdacp = [ Cp_array(1:B-2) Cp_array(B-1) 0 ]/sum(Cp_array)*vfs;
% vdacn = [ Cn_array(1:B-2)             0 0 ]/sum(Cn_array)*vfs;


%%%%%%%%%%%%
% Cp_array = [280.623 140.4002 70.2026 35.1031 17.5516 8.77588 4.38217 2.18867 1.10241 0.550779 0.550779];
% Cn_array = [280.623 140.4002 70.2026 35.1031 17.5516 8.77588 4.38217 2.18867 1.10241 0.550779 0.550779];
% Cu_mu = 0.550779;       % Cu = Cu_mu (fF)
% Cp_array = [ 2.^[(B-3):-1:0 0] ];
% Cn_array = [ 2.^[(B-3):-1:0 0] ];
% Cu_mu = 1;       % Cu = Cu_mu (fF)
% Cp_array = [70.4 35.2 17.6 8.8 4.4 2.2 1.1 0.55 0.55];
% Cn_array = [70.4 35.2 17.6 8.8 4.4 2.2 1.1 0.55 0.55];
% Cu_mu = 0.55;       % Cu = Cu_mu (fF)
% Cp_array = [70.2026 35.1031 17.5516 8.77588 4.38217 2.18867 1.10241 0.550779 0.550779];
% Cn_array = [70.2026 35.1031 17.5516 8.77588 4.38217 2.18867 1.10241 0.550779 0.550779];
% Cu_mu = 0.550779;       % Cu = Cu_mu (fF)
% Cp = 40;                % Parasitic Cap = Cp (fF)
% 
% Cu_sigma = 0.0000; % Cap mismatch variance relative to Cu
% C_nums = 2.^[(B-3):-1:0 0];
% for i = 1:length(C_nums)
%     Cp_array(1,i) = sum(randn(1,C_nums(i)))*Cu_sigma + Cp_array(1,i);
%     Cn_array(1,i) = sum(randn(1,C_nums(i)))*Cu_sigma + Cn_array(1,i);
% end
% 
% vdacp = [ Cp_array(1:B-2) Cp_array(B-1) 0 ]/(sum(Cp_array)+Cp)*vfs;
% vdacn = [ Cn_array(1:B-2)             0 0 ]/(sum(Cn_array)+Cp)*vfs;



% Cp_array = [ 2.^[(B-2):-1:0 0] ];
% Cn_array = [ 2.^[(B-2):-1:0 0] ];
% Cp_array = [140.2287 70.2026 35.1031 17.55156 8.77588 4.38217 2.18867 1.10241 0.550779 0.550779];
% Cn_array = [140.2287 70.2026 35.1031 17.55156 8.77588 4.38217 2.18867 1.10241 0.550779 0.550779];
% Cp_array = [140.2451 70.2098 35.1054 17.55156 8.77588 4.38217 2.18867 1.10241 0.550779 0.550779];
% Cn_array = [140.2451 70.2098 35.1054 17.55156 8.77588 4.38217 2.18867 1.10241 0.550779 0.550779];
% Cu_mu = 0.550779;       % Cu = Cu_mu (fF)
Cu_mu = 1;       % Cu = Cu_mu (fF)
Cp = 0;                % Parasitic Cap = Cp (fF)

Cu_sigma = 0.0000; % Cap mismatch variance relative to Cu
C_nums = 2.^[(B-2):-1:0 0];
for i = 1:length(C_nums)
    Cp_array(1,i) = sum(randn(1,C_nums(i)))*Cu_sigma + Cp_array(1,i);
    Cn_array(1,i) = sum(randn(1,C_nums(i)))*Cu_sigma + Cn_array(1,i);
end

vdacp = [ Cp_array(1:B-2) Cp_array(B-1) 0 ]/(sum(Cp_array)+Cp)*vfs;
vdacn = [ Cn_array(1:B-2) Cp_array(B-1) 0 ]/(sum(Cn_array)+Cp)*vfs;

%% sampling setting
fB=50e6; 
fnyq=2*fB; Tnyq=1/fnyq; %% nyquist output frequency
fclk = fnyq * B; Tclk = 1/fclk; 
N = 1000;

% input sine-signal setting
fsig = 17.3e6;
vin_amp_dB=-0.5; vin_amp = vfs*10^(vin_amp_dB/20); % single-ended amplitude

% simulation setting
t_sim = Tnyq*N;
options=simset('RelTol', 1e-3, 'MaxStep', 1/fclk);

% start simulation
open('SAR.slx'); 
out = sim('SAR', t_sim, options);

% extract result
sar_dout = out.SAR_DOUT(end-N+1:end,1)';
sar_dacp = out.SAR_DACP(end-N+1:end,1)';
sar_dacn = out.SAR_DACN(end-N+1:end,1)';


% calculate snr
fbin_sig=fsig/fnyq; w=hann_pv(N); fbin_L=3; fbin_H=N*fB/fnyq;
[snr, ptot]      = calcSNR(sar_dout,fbin_sig,fbin_L,fbin_H,w,N);


figure(1);
plot(sar_dout);

figure(2); %clf; 
plot(linspace(0,fnyq/2,N/2), ptot(1:N/2), 'r','linewidth',1.0); 
text_handle= text(floor(fB/8),-50, sprintf('SQNR = %4.1fdB',snr),'FontSize',18);
xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]'); title('Incremental ADC'); axis([0 fnyq/2 -130 0]);
title ([ 'SAR;  V_P= ',num2str(vin_amp_dB),'dBFS;  ','Freq=',num2str(fsig/1e3),'kHz;  ' ,'Conversion Rate=',num2str(fnyq/1e6),'MHz; ' ]);
legend('SQNR');


i_start_sar_dac = 12;
i_end_sar_dac = i_start_sar_dac+B-1;
figure(3) %dac voltage transition of Merged Capacitor Swiching (MCS) scheme
hold on;
stairs([sar_dacp(i_start_sar_dac:i_end_sar_dac) sar_dacp(i_end_sar_dac)], 'b');
stairs([sar_dacn(i_start_sar_dac:i_end_sar_dac) sar_dacn(i_end_sar_dac)], 'r');
hold off;
legend('DAC_n','DAC_p');
xlabel("(Nth Comparison Cycle)"); ylabel("(V)"); 


figure(4) %dac voltage transition of Monotonic Swiching (MS) scheme
hold on;
delta_voltage = vfs * 2.^(-[0:1:B-1]);
sar_dacp_MS = sar_dacp(i_start_sar_dac);
sar_dacn_MS = sar_dacn(i_start_sar_dac);
for i = 1:B-1
    if(sar_dacp_MS(i) > sar_dacn_MS(i))
        sar_dacp_MS(i+1) = sar_dacp_MS(i) - delta_voltage(i);
        sar_dacn_MS(i+1) = sar_dacn_MS(i);
    else
        sar_dacp_MS(i+1) = sar_dacp_MS(i);
        sar_dacn_MS(i+1) = sar_dacn_MS(i) - delta_voltage(i);
    end
end
stairs([sar_dacp_MS sar_dacp_MS(end)], 'b');
stairs([sar_dacn_MS sar_dacn_MS(end)], 'r');
hold off;
legend('DAC_n','DAC_p');
xlabel("(Nth Comparison Cycle)"); ylabel("(V)"); 


figure(5) % digital axis code-range = [0,1,2,3,4,5,6,7]/analog axis Voltage-range = [-1.1]
plot([-1:0.25:1],[2,2,2,2,2,2,2,2,2],'Color','black','Marker','+','MarkerEdgeColor','black');

figure(6) %3-Bit dac voltage transition of Merged Capacitor Swiching (MCS) scheme
hold on;
delta_voltage = vfs * 2.^(-[0:1:B-1]);
sar_dacp_MCS = 0.65;
sar_dacn_MCS = -0.65;
for i = 1:3
    if(sar_dacp_MCS(i) > sar_dacn_MCS(i))
        sar_dacp_MCS(i+1) = sar_dacp_MCS(i) - delta_voltage(i)/2;
        sar_dacn_MCS(i+1) = sar_dacn_MCS(i) + delta_voltage(i)/2;
    else
        sar_dacp_MCS(i+1) = sar_dacp_MCS(i) + delta_voltage(i)/2;
        sar_dacn_MCS(i+1) = sar_dacn_MCS(i) - delta_voltage(i)/2;
    end
end
stairs([sar_dacp_MCS sar_dacp_MCS(end)], 'b');
stairs([sar_dacn_MCS sar_dacn_MCS(end)], 'r');
hold off;
legend('DAC_n','DAC_p');
xlabel("(Nth Comparison Cycle)"); ylabel("(V)"); 

figure(7) %3-Bit dac voltage transition of Monotonic Swiching (MS) scheme
hold on;
delta_voltage = vfs * 2.^(-[0:1:B-1]);
sar_dacp_MS = 0.65;
sar_dacn_MS = -0.65;
for i = 1:3
    if(sar_dacp_MS(i) > sar_dacn_MS(i))
        sar_dacp_MS(i+1) = sar_dacp_MS(i) - delta_voltage(i);
        sar_dacn_MS(i+1) = sar_dacn_MS(i);
    else
        sar_dacp_MS(i+1) = sar_dacp_MS(i);
        sar_dacn_MS(i+1) = sar_dacn_MS(i) - delta_voltage(i);
    end
end
stairs([sar_dacp_MS sar_dacp_MS(end)], 'b');
stairs([sar_dacn_MS sar_dacn_MS(end)], 'r');
hold off;
legend('DAC_n','DAC_p');
xlabel("(Nth Comparison Cycle)"); ylabel("(V)"); 
