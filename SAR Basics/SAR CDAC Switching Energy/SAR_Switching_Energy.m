clear;
%% User Configuration
B_SAR = 10;


%% Parameters for Conventional Switching (CS)
Energy_per_switching_CS = [];
C_array_CS = 2.^([(B_SAR-1):-1:0 0]);
CA_CB_coefficients_CS = [];
Q_CA_CB_RedistributionTable_CS = [];
Average_SW_Energy_each_bit_CS = 0;
Average_SW_Energy_CS = 0;
Energy_of_each_Dcode_CS = [];
% stat_cs = zeros(1,B_SAR);
for i = 1:B_SAR+1
    CA = C_array_CS(i);
    CB = 2^(B_SAR)-C_array_CS(i);
    Ctot = 2^(B_SAR);
    CA_CB_coefficients_CS(1,i) = CA*CB/Ctot;
end  

for nth_bit = 1:B_SAR
    CA = C_array_CS(nth_bit);
    CB = 2^(B_SAR)-C_array_CS(nth_bit);
    Ctot = 2^(B_SAR);
    
    CB_array_CS = C_array_CS;
    CB_array_CS(nth_bit) = 0;
    Q_redistribution = -CA_CB_coefficients_CS(1,nth_bit).*[CB_array_CS./CB];
    Q_redistribution(nth_bit) = CA_CB_coefficients_CS(1,nth_bit);
    
    Q_CA_CB_RedistributionTable_CS(nth_bit,:) = Q_redistribution;
end
%% Switching Energy of Conventional Switching (CS)
DACP_CS_pre = zeros(1,B_SAR+1); 
DACN_CS_pre = ones(1,B_SAR+1);  
DACP_CS_cur = zeros(1,B_SAR+1); DACP_CS_cur(1) = 1;
DACN_CS_cur = ones(1,B_SAR+1);  DACN_CS_cur(1) = 0;
E_pre_switch = 2*CA_CB_coefficients_CS(1);
E_cur = E_pre_switch;

[Energy_per_switching_CS, Energy_of_each_Dcode_CS] = sar_CS_switching(1,B_SAR,1,DACP_CS_pre,DACN_CS_pre,DACP_CS_cur,DACN_CS_cur,[E_pre_switch],Q_CA_CB_RedistributionTable_CS,0,E_cur,[]); % (step,B_SAR,)

Average_SW_Energy_each_bit_CS = 0;
for i =0:B_SAR-1
    i_start = 2^i;
    i_end = i_start*2-1;
    temp = Energy_per_switching_CS(i_start:i_end);
    Average_SW_Energy_each_bit_CS(i+1) = mean(temp);
end
 Average_SW_Energy_CS= sum(Average_SW_Energy_each_bit_CS);


 

%% Parameters for Monotonic Switching (MS)
Energy_per_switching_MS = [];
C_array_MS = 2.^([(B_SAR-2):-1:0 0]);
CA_CB_coefficients_MS = [];
Q_CA_CB_RedistributionTable_MS = [];
Average_SW_Energy_each_bit_MS = 0;
Average_SW_Energy_MS = 0;
Energy_of_each_Dcode_MS = [];
% stat_cs = zeros(1,B_SAR);
for i = 1:B_SAR
    CA = C_array_MS(i);
    CB = 2^(B_SAR-1)-C_array_MS(i);
    Ctot = 2^(B_SAR-1);
    CA_CB_coefficients_MS(1,i) = CA*CB/Ctot;
end  

for nth_bit = 1:(B_SAR-1)
    CA = C_array_MS(nth_bit);
    CB = 2^(B_SAR-1)-C_array_MS(nth_bit);
    Ctot = 2^(B_SAR-1);
    
    CB_array_MS = C_array_MS;
    CB_array_MS(nth_bit) = 0;
    Q_redistribution = -CA_CB_coefficients_MS(1,nth_bit).*[CB_array_MS./CB];
    Q_redistribution(nth_bit) = CA_CB_coefficients_MS(1,nth_bit);
    
    Q_CA_CB_RedistributionTable_MS(nth_bit,:) = Q_redistribution;
end
%% Switching Energy of Monotonic Switching (MS)
DACP_MS_pre = ones(1,B_SAR); 
DACN_MS_pre = ones(1,B_SAR);  
DACP_MS_cur = DACP_MS_pre;
DACN_MS_cur = DACN_MS_pre;
E_pre_switch = 0;
E_cur = E_pre_switch;

[Energy_per_switching_MS, Energy_of_each_Dcode_MS] = sar_MS_switching(1,B_SAR,1,DACP_MS_pre,DACN_MS_pre,DACP_MS_cur,DACN_MS_cur,[E_pre_switch],Q_CA_CB_RedistributionTable_MS,0,E_cur,[]); % (step,B_SAR,)

Average_SW_Energy_each_bit_MS = 0;
for i =0:B_SAR-1
    i_start = 2^i;
    i_end = i_start*2-1;
    temp = Energy_per_switching_MS(i_start:i_end);
    Average_SW_Energy_each_bit_MS(i+1) = mean(temp);
end
 Average_SW_Energy_MS= sum(Average_SW_Energy_each_bit_MS);


 

%% Parameters for Merged Capacitor Switching (MCS)
Energy_per_switching_MCS = [];
C_array_MCS = 2.^([(B_SAR-2):-1:0 0]);
CA_CB_coefficients_MCS = [];
Q_CA_CB_RedistributionTable_MCS = [];
Average_SW_Energy_each_bit_MCS = 0;
Average_SW_Energy_MCS = 0;
Energy_of_each_Dcode_MCS = [];
% stat_cs = zeros(1,B_SAR);
for i = 1:B_SAR
    CA = C_array_MCS(i);
    CB = 2^(B_SAR-1)-C_array_MCS(i);
    Ctot = 2^(B_SAR-1);
    CA_CB_coefficients_MCS(1,i) = CA*CB/Ctot;
end  

for nth_bit = 1:(B_SAR-1)
    CA = C_array_MCS(nth_bit);
    CB = 2^(B_SAR-1)-C_array_MCS(nth_bit);
    Ctot = 2^(B_SAR-1);
    
    CB_array_MCS = C_array_MCS;
    CB_array_MCS(nth_bit) = 0;
    Q_redistribution = -CA_CB_coefficients_MCS(1,nth_bit).*[CB_array_MCS./CB];
    Q_redistribution(nth_bit) = CA_CB_coefficients_MCS(1,nth_bit);
    
    Q_CA_CB_RedistributionTable_MCS(nth_bit,:) = Q_redistribution;
end
%% Switching Energy of Merged Capacitor Switching (MCS)
DACP_MCS_Vcm_pre = 0.5*ones(1,B_SAR); 
DACN_MCS_Vcm_pre = 0.5*ones(1,B_SAR);  
DACP_MCS_Vref_pre = zeros(1,B_SAR); 
DACN_MCS_Vref_pre = zeros(1,B_SAR); 
DACP_MCS_Vcm_cur = DACP_MCS_Vcm_pre;
DACN_MCS_Vcm_cur = DACN_MCS_Vcm_pre;
DACP_MCS_Vref_cur = DACP_MCS_Vref_pre;
DACN_MCS_Vref_cur = DACN_MCS_Vref_pre; 

E_pre_switch = 0;
E_cur = E_pre_switch;

[Energy_per_switching_MCS, Energy_of_each_Dcode_MCS] = sar_MCS_switching(1,B_SAR,1,DACP_MCS_Vcm_pre,DACN_MCS_Vcm_pre,DACP_MCS_Vcm_cur,DACN_MCS_Vcm_cur,DACP_MCS_Vref_pre,DACN_MCS_Vref_pre,DACP_MCS_Vref_cur,DACN_MCS_Vref_cur,[E_pre_switch],Q_CA_CB_RedistributionTable_MCS,0,E_cur,[]); % (step,B_SAR,)

Average_SW_Energy_each_bit_MCS = 0;
for i =0:B_SAR-1
    i_start = 2^i;
    i_end = i_start*2-1;
    temp = Energy_per_switching_MCS(i_start:i_end);
    Average_SW_Energy_each_bit_MCS(i+1) = mean(temp);
end
 Average_SW_Energy_MCS= sum(Average_SW_Energy_each_bit_MCS);
 
 
 figure(1);
set(gcf, 'color', [1 1 1]);
 hold on;
 plot(Energy_of_each_Dcode_CS,'LineWidth',2);
 plot(Energy_of_each_Dcode_MS,'LineWidth',2);
 plot(Energy_of_each_Dcode_MCS,'LineWidth',2);
 hold off;
 xlim([1,2^B_SAR]);
set(gca, 'FontSize', 12);
xlabel('Code(s)','FontSize',14);
ylabel('Switching Energy (CV_r_e_f ^2)','FontSize',14);
%  yticks([100 400 800 1200 1600 2000]);
legend({'CS','MS','MCS'},'FontSize',16);

Normalized_Eavg_CS  = Average_SW_Energy_CS /Average_SW_Energy_CS;
Normalized_Eavg_MS  = Average_SW_Energy_MS /Average_SW_Energy_CS;
Normalized_Eavg_MCS = Average_SW_Energy_MCS/Average_SW_Energy_CS;
