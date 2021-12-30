%% Incremental ADC-  AC test
clear
addpath ('D:\LocalLaptop\õìëË\å§ãÜèä\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\õìëË\å§ãÜèä\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');

%% User Parameters
opgain_dB_ARR_USER = 80;
% opgain_dB_ARR_USER = [70:10:100];
% SAR_GAIN_ERROR_ARR_USER = [0.94:0.01:1.06];
SAR_GAIN_ERROR_ARR_USER = [0.98:0.005:1.1];
% vin_amp_dB_ARR_USER = [-80:5:-10, -9.5:0.5:0];  
vin_amp_dB_ARR_USER = [-3];  
Compen_Gain_Error_ARR_USER=-20:100;

%% Modulator Parameters
%vfs=3.3; vref=vfs/2; vcm=vfs/2; vicm=0; %vicm=vcm;  % full-scale vfs=3.3V
vref=0.9; vfs=vref*2; VFS=vfs; %vcm=vfs/2; vicm=vcm; % full-scale vfs=2V
vin_amp_dB_ARR=vin_amp_dB_ARR_USER;
vin_amp=10.^(vin_amp_dB_ARR/20)*vref;

%% SAR-PART Parameters
SAR_B=10;

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
opgain_dB_ARR=opgain_dB_ARR_USER; opgain=10.^(opgain_dB_ARR/20); 

%% ExCount Gain Error Cal. Parameters
Compen_Gain_Error_ARR=Compen_Gain_Error_ARR_USER;


%% Open Simulink diagram
options=simset('RelTol', 1e-3, 'MaxStep', 1/fs);
%% model selection for different Simulink versions
open('ExCount_wi_SarGainError.mdl');

ExCount_arr1 = []; 
IADC_pp_out_arr1 = []; 
SAR_out_arr1 = []; 
for i=1:length(opgain)
    for j=1:length(SAR_GAIN_ERROR_ARR)
        A = opgain(i);
        SAR_GAIN_ERROR = SAR_GAIN_ERROR_ARR(j);
    %     alfa = alfa_arr(i);
        sim('ExCount_wi_SarGainError', t_sim, options); % 2015 Simulink
        ExCount_arr1(i,j,:,:) = out_no_cal(end-n_nyq+1:end,:)';
        IADC_pp_out_arr1(i,j,:,:) = IADC_pp_out(end-n_nyq+1:end,:)'; 
        SAR_out_arr1(i,j,:,:) = SAR_out(end-n_nyq+1:end,:)'; 
    end
end

%% FFT @ Nyquist Rate
sqnr_ExCount_out = zeros(length(opgain),length(SAR_GAIN_ERROR_ARR),length(vin_amp_dB_ARR));
ptot_ExCount_out = zeros(length(opgain),length(SAR_GAIN_ERROR_ARR),length(vin_amp_dB_ARR),n_nyq);

sqnr_ExCount_out_cal_by_div = zeros(length(opgain),length(SAR_GAIN_ERROR_ARR),length(vin_amp_dB_ARR),length(Compen_Gain_Error_ARR));
ptot_ExCount_out_cal_by_div  = zeros(length(opgain),length(SAR_GAIN_ERROR_ARR),length(vin_amp_dB_ARR),length(Compen_Gain_Error_ARR),n_nyq);
sqnr_ExCount_out_cal_by_subtract = zeros(length(opgain),length(SAR_GAIN_ERROR_ARR),length(vin_amp_dB_ARR),length(Compen_Gain_Error_ARR));
ptot_ExCount_out_cal_by_subtract  = zeros(length(opgain),length(SAR_GAIN_ERROR_ARR),length(vin_amp_dB_ARR),length(Compen_Gain_Error_ARR),n_nyq);


fbin_sig=fsig/fnyq; w=hann_pv(n_nyq); fbin_L=3; fbin_H=n_nyq*fB/fnyq;
for i = 1:length(opgain)
    for j = 1:length(SAR_GAIN_ERROR_ARR)
        for k=1:length(vin_amp_dB_ARR)
%             ExCount_out1= reshape(ExCount_arr1(i,j,k,:),[1,n_nyq]);
%             [sqnr_ExCount_out(i,j,k),ptot_ExCount_out(i,j,k,:)]   = calcSNR(ExCount_out1,fbin_sig,fbin_L,fbin_H,w,n_nyq);
                
            
            IADC_out = reshape(IADC_pp_out_arr1(i,j,k,:),[1,n_nyq]);
            SAR_out = reshape(SAR_out_arr1(i,j,k,:),[1,n_nyq]);
            %%** A. W/O SAR-gain-error calibration
            ExCount_out2 = IADC_out+SAR_out; 
            [sqnr_ExCount_out(i,j,k),ptot_ExCount_out(i,j,k,:)]   = calcSNR(ExCount_out2,fbin_sig,fbin_L,fbin_H,w,n_nyq);
      
               
            %%** B. W/ SAR-gain-error calibration
            %     %%**** B-1 (Devision)
%             sqnr_test3_array=[];
            for i_gain_error_cal = 1:length(Compen_Gain_Error_ARR)
                SAR_res3 = (SAR_out*2^(SAR_B) / (2^(SAR_B)+Compen_Gain_Error_ARR(i_gain_error_cal)));
                ExCount_out3 =  IADC_out + SAR_res3;
                ExCount_out3 = ExCount_out3*(2/(M-1)/(M-2));

                [sqnr_test3,ptot_test3] = calcSNR(ExCount_out3, fbin_sig, fbin_L, fbin_H, w, n_nyq);
                sqnr_ExCount_out_cal_by_div(i,j,k,i_gain_error_cal) = sqnr_test3;
                ptot_ExCount_out_cal_by_div(i,j,k,i_gain_error_cal,:) = ptot_test3;
                
%                 [sqnr_test3, ptot_test3] = calcSNR(ExCount_out3, fbin_sig, fbin_L, fbin_H, w, N);
%                 A = sqnr_test3;
%                 sqnr_test3_array(1,end+1) = sqnr_test3;
            end
            %%**** B-2 (Subtraction)
%             sqnr_test4_array=[];
            for i_gain_error_cal = 1:length(Compen_Gain_Error_ARR)
                temp_sub = (SAR_out*2^(SAR_B))*Compen_Gain_Error_ARR(i_gain_error_cal) / (2^(SAR_B));
                SAR_res4 = ( ( (SAR_out*2^(SAR_B)) -temp_sub) / (2^(SAR_B)));
                ExCount_out4 = IADC_out+SAR_res4;   
                ExCount_out4 = ExCount_out4*(2/(M-1)/(M-2));

                [sqnr_ExCount_out_cal_by_subtraction(i,j,k,i_gain_error_cal),ptot_ExCount_out_cal_by_subtraction(i,j,k,i_gain_error_cal,:)] = calcSNR(ExCount_out4, fbin_sig, fbin_L, fbin_H, w, n_nyq);
%                 [sqnr_test4, ptot_test4] = calcSNR(yt, fbin_sig, fbin_L, fbin_H, w, N);
%                 A = sqnr_test4
%                 sqnr_test4_array(1,end+1) = sqnr_test4;
            end
        end
        
        

%     
%     %%**** B-2 (Subtraction)
%     sqnr_test4_array=[];
%     for i_gain_error_cal = 1:length(Compen_Gain_Error_ARR)
%         temp_sub = D2_Array(i_results,:) *Compen_Gain_Error_ARR(i_gain_error_cal) / (2^(SAR_B));
%         sar_res = ( (D2_Array(i_results,:) -temp_sub) / (2^(SAR_B)));
%         yt = D1_Array(i_results,:);
%         yt = yt+sar_res(1:n_nyq);
%         % yt = yt+residual(1:n_nyq)
%         yt = yt*(2/M/(M-1));
% 
%         fbin_sig=fsig/fnyq; w=hann_pv(n_nyq); fbin_L=3; fbin_H=n_nyq*fB/fnyq; N = n_nyq;
%         [sqnr_test4, ptot_test4] = calcSNR(yt, fbin_sig, fbin_L, fbin_H, w, N);
%         A = sqnr_test4
%         sqnr_test4_array(1,end+1) = sqnr_test4;
%     end
% 
%     figure(1+i_results); clf;
%     hold on;
%     scatter(Compen_Gain_Error_ARR,sqnr_test3_array,[],'green');
%     scatter(Compen_Gain_Error_ARR,sqnr_test4_array,[],'magenta');
%     peak_sqnr_test3 = max(sqnr_test3_array);
%     peak_sqnr_test4 = max(sqnr_test4_array);
%     i_peak_sqnr_test3 = find(sqnr_test3_array==peak_sqnr_test3);
%     i_peak_sqnr_test4 = find(sqnr_test4_array==peak_sqnr_test4);
%     value_compen_test3 = Compen_Gain_Error_ARR(i_peak_sqnr_test3);
%     value_compen_test4 = Compen_Gain_Error_ARR(i_peak_sqnr_test4);
%     % text_handle= text(11,min(sqnr_test3_array)+2, sprintf('Peak SQNR = %4.1fdB (x=%2d)',peak_sqnr_test3,value_compen_test3),'Color','green','FontSize',14,'HorizontalAlignment','right');
%     % text_handle= text(11,min(sqnr_test3_array)+1.5, sprintf('Peak SQNR = %4.1fdB (x=%2d)',peak_sqnr_test4,value_compen_test4),'Color','magenta','FontSize',14,'HorizontalAlignment','right');
%     text_handle= text(mean(Compen_Gain_Error_ARR),min(sqnr_test3_array)+2, sprintf('Peak SQNR = %4.1fdB (x=%2d)',peak_sqnr_test3,value_compen_test3),'Color','green','FontSize',14,'HorizontalAlignment','center');
%     text_handle= text(mean(Compen_Gain_Error_ARR),min(sqnr_test3_array)+1, sprintf('Peak SQNR = %4.1fdB (x=%2d)',peak_sqnr_test4,value_compen_test4),'Color','magenta','FontSize',14,'HorizontalAlignment','center');
%     legend('Division','Subtraction','Location','northwest');
%     xlim([Compen_Gain_Error_ARR(1),Compen_Gain_Error_ARR(end)]);
%     xlabel('x'); ylabel('SQNR [dB]');  
%     title ({ ['Extended Counting',' (',num2str(SAR_B),'-Bit SAR, ',num2str(SAR_N_Comparison),' Quantization Cycles)']; ['V_P = ',num2str(sig_dB),' dBFS;  ','Fsig = ',num2str(fsig/1e3),' kHz;  ' ,'Fs = ',num2str(fs/1e6),' MHz; ','OSR = ', num2str(OSR)] });
% 
%     i_peak_sqnr_Devision_Array    = [i_peak_sqnr_Devision_Array    i_peak_sqnr_test3];
%     i_peak_sqnr_Subtraction_Array = [i_peak_sqnr_Subtraction_Array i_peak_sqnr_test4];
%     
    
    
    
    end
end

%% SNR v.s. Amplitude (Incremental ADC)
line_type=["--r*","--y*","--g*","--b*","--m*"];

figure (1); clf; %%
set(gcf, 'color', [1 1 1]);
hold on;
% lengend_txt={};
for i_opgain = 1:length(opgain)
    out_size = size(sqnr_ExCount_out);
    sqnr_ExCount_out_for_specific_opgain =  reshape(sqnr_ExCount_out(i_opgain,:,:),[length(SAR_GAIN_ERROR_ARR),length(vin_amp)]);
    plot(SAR_GAIN_ERROR_ARR,max(sqnr_ExCount_out_for_specific_opgain,[],2)','linewidth',0.5,'DisplayName',sprintf("A_o = %4.1fdB",opgain_dB_ARR(i_opgain)) ); 
%     txt = sprintf("A_o=%4.1fdB",opgain_dB_ARR);
%     lengend_txt(end+1) = {txt};
end
% legend(lengend_txt);
legend('FontSize',14);
legend('boxoff');
title ({['(W/O Calibration)'],['peak SQNR among all imput amplitude V.S.  SAR\_GAIN\_ERROR \alpha']});
xlabel ('SAR\_GAIN\_ERROR \alpha','FontSize',32); ylabel ('SQNR (dB)','FontSize',32); 
hold off;
set(gca, 'FontSize',11);



opgain_dB_chosed = 80;
Vsig_dB_chosed = -3;
i_opgain_dB_ARR_chosed = find(opgain_dB_ARR==opgain_dB_chosed);
i_vin_amp_dB_ARR_chosed  = find(vin_amp_dB_ARR==Vsig_dB_chosed);
figure (2); clf;
set(gcf, 'color', [1 1 1]);
hold on;
for i_SAR_Gain_ERROR = 1:length(SAR_GAIN_ERROR_ARR)
        sqnr_ExCount_out_NO_cal =  reshape(sqnr_ExCount_out(i_opgain,:,i_vin_amp_dB_ARR_chosed),[1,length(SAR_GAIN_ERROR_ARR)]);
end
scatter(SAR_GAIN_ERROR_ARR,sqnr_ExCount_out_NO_cal,'blue');
title ({['(W/O Calibration) SQNR  V.S.  SAR\_GAIN\_ERROR \alpha'];['(@V_P= ',num2str(Vsig_dB_chosed),'dBFS; OP GAIN=',num2str(opgain_dB_chosed),'dB) ']});
xlabel ('SAR\_GAIN\_ERROR \alpha','FontSize',32); ylabel ('SQNR (dB)','FontSize',32); 
hold off;
set(gca, 'FontSize',12);


figure (3); clf;
set(gcf, 'color', [1 1 1]);
peak_sqnr_cal_by_div = [];
peak_sqnr_cal_by_subtraction = [];
for i_SAR_Gain_ERROR = 1:length(SAR_GAIN_ERROR_ARR)
    tmp = sqnr_ExCount_out_cal_by_div(i_opgain_dB_ARR_chosed,i_SAR_Gain_ERROR,i_vin_amp_dB_ARR_chosed,:);
    sqnr_cal_by_div_SWEEP =reshape(tmp,[1,length(Compen_Gain_Error_ARR)]);
    peak_sqnr_cal_by_div(1,end+1) = max(sqnr_cal_by_div_SWEEP);
    
    tmp = sqnr_ExCount_out_cal_by_subtraction(i_opgain_dB_ARR_chosed,i_SAR_Gain_ERROR,i_vin_amp_dB_ARR_chosed,:);
    sqnr_cal_by_subtraction_SWEEP =reshape(tmp,[1,length(Compen_Gain_Error_ARR)]);
    peak_sqnr_cal_by_subtraction(1,end+1) = max(sqnr_cal_by_subtraction_SWEEP);
end
hold on;
scatter(SAR_GAIN_ERROR_ARR,peak_sqnr_cal_by_div,'green');
scatter(SAR_GAIN_ERROR_ARR,peak_sqnr_cal_by_subtraction,'magenta');
legend('By Division','By Subtraction','Location','northwest','FontSize',14);
legend('boxoff');
title ({['(W/ Calibration) SQNR  V.S.  SAR\_GAIN\_ERROR \alpha'];['(@V_P= ',num2str(Vsig_dB_chosed),'dBFS; OP GAIN=',num2str(opgain_dB_chosed),'dB) ']});
xlabel ('SAR\_GAIN\_ERROR \alpha','FontSize',32); ylabel ('SQNR (dB)','FontSize',32); 
hold off;
set(gca, 'FontSize',13);

%% END