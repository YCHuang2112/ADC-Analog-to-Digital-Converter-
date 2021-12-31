%% Incremental ADC-  AC test
clear
addpath ('D:\LocalLaptop\õìëË\å§ãÜèä\adc\Matlab_Toolbox\MIT Design Tools/HSPICE_Toolbox/HspiceToolbox');
addpath ('D:\LocalLaptop\õìëË\å§ãÜèä\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\õìëË\å§ãÜèä\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');
%% 9 Corner
%           Temp    PWRA    Corner
% Corner1   50      1.8     TT
% Corner2   0       1.62    SS
% Corner3   100     1.62    SS
% Corner4   0       1.98    SS
% Corner5   100     1.98    SS
% Corner6   0       1.62    FF
% Corner7   100     1.62    FF
% Corner8   0       1.98    FF
% Corner9   100     1.98    FF

Corner_Temp=[50 0 100 0 100 0 100 0 100];
Corner_PWRA=[1.8 1.62 1.62 1.98 1.98 1.62 1.62 1.98 1.98];
Corner_Name={'TT' 'SS' 'SS' 'SS' 'SS' 'FF' 'FF' 'FF' 'FF'};
Corner_color = [[13  105 255];[0   153 255];[0   180 200];[72  234 195];[0   255 204];
           [1   255 98 ];[173 212 6  ];[237 205 47 ];[86  171 29 ];[0   200 0  ];
           [213 120 51 ];[250 10  10 ];[242 38  208];[167 75  251];[100 0   200];
           [180 0   180];]/255;
% Corner_color=[];
% for i=1:16
%     Corner_color=[Corner_color;[4*i   15*i   255-15*i]/255;]; %Blue
% %     Corner_color=[Corner_color;[245 245 245]/255;]; %White Smoke
% %     Corner_color=[Corner_color;[211 211 211]/255;]; %Light Gray
% %     Corner_color=[Corner_color;[192 192 192]/255;]; %Silver
% %     Corner_color=[Corner_color;[169 169 169]/255;]; %Dark Gray
% end
% Corner_color = winter(16);

%********************************************************************************************
%********************************* START: Parameter settings ********************************
%% User Parameters
%***** A-1. Group S1
N_points_for_analysis = 80; fB_user=1e6;   fsig_user=325e3; OSR_user=16; RST_PATTERN=1;

%***** B. (lowV,nomV,highV)=(1,62,1.8,1.98)
vref_9corner_array=Corner_PWRA/2;

%***** C.
User_Gain_Error_Compensation_Mode=2; % {1:devision; 2:subtraction; default:No-Compensation}

%% Modulator Parameters
vref=vref_9corner_array; vfs=vref*2; VFS=vfs; 
OSR=OSR_user; M=OSR;

%% ExCount Gain Error Cal. Parameters
Compen_Gain_Error_Array=0:60;
Division_1__Subtraction_2=User_Gain_Error_Compensation_Mode; % {1:devision; 2:subtraction; default:No-Compensation}

%% SAR-PART Parameters
SAR_B=10;

%% Input Signal parameters
fsig=fsig_user; Tsig=1/fsig; Vsig_dB = -3;


%% Frequency of Clock
fB=fB_user; fnyq=2*fB; T_nyq=1/fnyq; %% nyquist output frequency
fs=OSR*fnyq; Ts=1/fs;  
fclk=fs*1; tclk=Ts;

%% Numbers of nyquist points
n_nyq=N_points_for_analysis;
%*********************************  END: Parameter settings  ********************************
%********************************************************************************************

%********************************************************************************************
%********************************** START: Readin-Files Settings ****************************
%% Set Input Signal
filename_extensions = {'.tr0','.tr1','.tr2','.tr3','.tr4','.tr5','.tr6','.tr7','.tr8'};


%*** Group S1 (default:10Comp)
%****** S1-9Comp
SAR_N_Comparison = 9; 
sig_cell1 = {'v_sar_rst_b'};

%%% Option 1 (SAR_SAMPLING_Time=1_subperiod)
RST_PATTERN=3; filename_prefix = 'TOP_ADC_ExCount_v5/f2_cp2_v5';  process_code='Presim'; Compen_Gain_Error=38; Internal_Gain=2; sig_dB=-3; fsig=725e3;
% RST_PATTERN=3; filename_prefix = 'TOP_ADC_ExCount_v8/f2_cp2_1_period';  process_code='Presim'; Compen_Gain_Error=38; Internal_Gain=2; sig_dB=-3; fsig=725e3;

%%% Option 2 (SAR_SAMPLING_Time=1.5_subperiod)
% RST_PATTERN=1; filename_prefix = 'TOP_ADC_ExCount_v6/f2_cp2_v6';  process_code='Presim'; Compen_Gain_Error=38; Internal_Gain=2; sig_dB=-3; fsig=725e3;
% RST_PATTERN=1; filename_prefix = 'TOP_ADC_ExCount_v8/f2_cp2_1.5_period';  process_code='Presim'; Compen_Gain_Error=38; Internal_Gain=2; sig_dB=-3; fsig=725e3;


sig_cell3 = {'v_sa<0>','v_sa<1>','v_sa<2>','v_sa<3>','v_sa<4>','v_sa<5>','v_sa<6>','v_sa<7>','v_sa<8>'};


%*********************************  END: Readin-Files Settings  *****************************
%********************************************************************************************

%********************************************************************************************
%********************************* START: Post-Processing Loop  *****************************
sig_cell2 = {'v_ds_out<0>','v_ds_out<1>','v_ds_out<2>','v_ds_out<3>'};

sqnr_test_results=[];ptot_test_results=[];
D1_Array=[];D2_Array=[];Residue_Array=[];
for i_test = 1:length(filename_extensions)
    %********************************************************************************************
    %********************************** START: Signals Readin ***********************************
    %% Read Signals
    %***** A. load signals from a file
    x = loadsig([filename_prefix,filename_extensions{i_test}]);
    sig = lssig(x);

    %***** B. load signal: reset
    rst_b = evalsig(x,sig_cell1{1});



    i_start = 1;
    for i=1:length(rst_b)
        if(rst_b(i) > vfs/2)
               i_start = i;
        else
                break;
        end
    end

    rst(1:16,1) = 0;
%     rst(2,1) = 0;
    if RST_PATTERN==1
        for i=16:length(rst_b)
            if( rst_b(i)<vref(i_test) && rst_b(i-1)<vref(i_test))
                rst(i,1)=1;
            else
                rst(i,1)=0;
            end
        end
    elseif RST_PATTERN==2
        for i=16:length(rst_b)
            if( rst_b(i)<vref(i_test) && rst_b(i-1)>vref(i_test))
                rst(i,1)=1;
            else
                rst(i,1)=0;
            end
        end
    elseif RST_PATTERN==3
        for i=16:length(rst_b)
            if( rst_b(i)>vref(i_test) && rst_b(i-1)<vref(i_test))
                rst(i,1)=1;
            else
                rst(i,1)=0;
            end
        end
    end
    for i=1:length(rst)
        if( rst(i)==1 )
            i_start=i+1;
            break;
        end
    end
    rst = rst(i_start:end);

    %***** C. load signal: IADC_douts
    clear y;
    for i = 1:length(sig_cell2)
        y(i,:) = evalsig(x,sig_cell2{i})';
    end
    for i=1:size(y,2)
        for j=1:size(y,1)
            if y(j,i)>0.9
                y(j,i)=1;
            else
                y(j,i)=-1;
            end
        end
    end
    test = sum(y,1)/4;
    test_pre = test;
    test = test(i_start:end)';

    %***** D. load signal: SAR_douts
    clear residue;
    for i = 1:length(sig_cell3)
        residue(i,:) = evalsig(x,sig_cell3{i})';
    end
    for i=1:size(residue,2)
        for j=1:size(residue,1)
            if residue(j,i)>0.9
                residue(j,i)=1;
            else
                residue(j,i)=-1;
            end
        end
    end
    sar_Quantization_cycles = length(sig_cell3);
    residue_len = size(residue,2);
    sar_res_pre = zeros(1,residue_len);
    for i = 1:sar_Quantization_cycles
        sar_res_pre = residue(i,:)*(2^(i-1)) + sar_res_pre;
    end

    sar_res_pre = sar_res_pre(i_start:end); 
    sar_res = [];
    if RST_PATTERN==1
        for i = 1:(length(sar_res_pre)-1)
            if rst(i) == 1
                sar_res = cat( 2,sar_res, sar_res_pre(i+1) );
            end 
        end
    elseif RST_PATTERN==2
        for i = 1:length(sar_res_pre)
            if rst(i) == 1
                sar_res = cat( 2,sar_res, sar_res_pre(i) );
            end 
        end
    elseif RST_PATTERN==3
        for i = 1:length(sar_res_pre)
            if rst(i) == 1
                sar_res = cat( 2,sar_res, sar_res_pre(i) );
            end 
        end
    end

    sar_res = sar_res(2:end);
    sar_res_save02 = sar_res;
    if(Division_1__Subtraction_2==1)
        %%**** Division
        sar_res = (sar_res / (2^(SAR_B)+Compen_Gain_Error));
    elseif(Division_1__Subtraction_2==2)
        %%**** Subtraction
        temp_sub = sar_res*Compen_Gain_Error / (2^(SAR_B));
        sar_res = ( (sar_res-temp_sub) / (2^(SAR_B)));
    else
        %%**** Default: No gain error compensation
        sar_res = (sar_res / (2^(SAR_B)));
    end


    %***** E. load signal: IADC Part's Residual
    residue_pre1 = (evalsig(x,'v_residualxg1o1_p')-evalsig(x,'v_residualxg1o1_n'))/VFS(i_test);
    residue_pre = residue_pre1(i_start:end);
    k = 1;
    for i=1:length(residue_pre)
        if(rst(i)==1)
            residue(k,1) = residue_pre(i);
            k = k+1;
        end
    end


    residue = residue * 1  ;


    %**********************************  END: Signals Readin  ***********************************
    %********************************************************************************************

    %********************************************************************************************
    %********************************** START: Post-Processing **********************************
    %% Set Simulink time infos
    t_sim = length(rst)/fclk;
    t_index=[1:length(rst)]'; t=(t_index-0.5)/fclk;

    %% Open Simulink diagram
    options=simset('RelTol', 1e-3, 'MaxStep', 1/fs);

    %% model selection for different Simulink versions
    open('iadc2_post_processing.mdl'); sim('iadc2_post_processing', t_sim, options); % 2015 Simulink

    %% Signal for Debug
%     XXX = [simulink_rst,[rst;0],simulink_input,[test;0],count1,count2];
    %**********************************  END: Post-Processing  **********************************
    %********************************************************************************************

    %********************************************************************************************
    %********************************* START: Spectrum Analysis *********************************
    %% DFT Analysis
    fbin_sig=fsig/fnyq; w=hann_pv(n_nyq); fbin_L=3; fbin_H=n_nyq*fB/fnyq; N = n_nyq;
    yt = zeros(1,n_nyq);
    for i = 1:n_nyq
        yt(i) =  idc2_pre(OSR*i,1)*1;
    end
    yt_pre = yt;
    yt = yt+residue(1:n_nyq);
    yt = yt*(2/M/(M-1));
    yt_fig1 = yt;

%     xxxx = [yt_pre',sar_res(1:N_points_for_analysis)'];
% 
%     XXX_sar = sar_res_save02(1:n_nyq) / (2^SAR_B);
%     XXX_IADC_res = residue(1:n_nyq);
%     ZZZ = XXX_sar./XXX_IADC_res;

    [sqnr_test, ptot_test] = calcSNR(yt, fbin_sig, fbin_L, fbin_H, w, N);
    
    sqnr_test_results  = [sqnr_test_results  sqnr_test];
    ptot_test_results = [ptot_test_results ptot_test];
    D1_Array=[D1_Array; yt_pre;];
    D2_Array=[D2_Array; sar_res_save02(1:n_nyq);];
    Residue_Array=[Residue_Array; residue(1:n_nyq);];
    %*********************************  END: Spectrum Analysis  *********************************
    %********************************************************************************************
end
%*********************************  END: Post-Processing Loop  ******************************
%********************************************************************************************

%********************************************************************************************
%********************************* START: Plotting Spectrum Analysis Results  ***************
%% A. DFT Analysis of (ExCount Dout + Unquantized Analog Residue)
figure(1); clf;
hold on;
lengend_txt={};
for i_results=1:size(ptot_test_results,2)
    plot(linspace(0,fnyq/2,n_nyq/2+1), ptot_test_results(1:n_nyq/2+1,i_results),'Color',Corner_color(i_results,:),'linewidth',2.0);
    text_handle= text(5,-1-i_results*5, sprintf('SQNR = %4.1f dB',sqnr_test_results(i_results)),'Color',Corner_color(i_results,:),'FontSize',12);
    txt = [sprintf('%5s',num2str(Corner_PWRA(i_results))),'V ',Corner_Name{i_results},sprintf(' %3s',num2str(Corner_Temp(i_results))),char(0176),'C',];
    lengend_txt(end+1) = {txt};
end
xlabel('Frequency [Hz]'); ylabel('PSD [dB/Hz]');  axis([0 fB -130 0]); 
title ({['Extended Counting with Unquantized Analog Residue (',process_code,')'];['V_P= ',num2str(Vsig_dB),'dBFS;  ','Fsig=',num2str(fsig/1e3),'kHz;  ' ,'Fs=',num2str(fs/1e6/OSR),'MHz; ','OSR= ', num2str(OSR) ]},'FontSize',14);

legend(lengend_txt);
hold off;

%% B. DFT Analysis of (ExCount Dout + Sweeping Compen_Gain_Error)
i_peak_sqnr_Devision_Array=[];i_peak_sqnr_Subtraction_Array=[];
for i_results=1:size(filename_extensions,2)
    %%**** B-1 (Devision)
    sqnr_test3_array=[];
    for i_gain_error_cal = 1:length(Compen_Gain_Error_Array)
    %     sar_res = (D2_Array(i_results,:)  / (2^(sar_bits)+Compen_Gain_Error_Array(i_gain_error_cal))) / Internal_Gain;
        sar_res = (D2_Array(i_results,:) / (2^(SAR_B)+Compen_Gain_Error_Array(i_gain_error_cal)));
        yt = D1_Array(i_results,:);
        yt = yt+sar_res(1:n_nyq);
    %     yt = yt+residual(1:n_nyq)
        yt = yt*(2/M/(M-1));

        fbin_sig=fsig/fnyq; w=hann_pv(n_nyq); fbin_L=3; fbin_H=n_nyq*fB/fnyq; N = n_nyq;
        [sqnr_test3, ptot_test3] = calcSNR(yt, fbin_sig, fbin_L, fbin_H, w, N);
        A = sqnr_test3
        sqnr_test3_array(1,end+1) = sqnr_test3;
    end
    
    %%**** B-2 (Subtraction)
    sqnr_test4_array=[];
    for i_gain_error_cal = 1:length(Compen_Gain_Error_Array)
        temp_sub = D2_Array(i_results,:) *Compen_Gain_Error_Array(i_gain_error_cal) / (2^(SAR_B));
        sar_res = ( (D2_Array(i_results,:) -temp_sub) / (2^(SAR_B)));
        yt = D1_Array(i_results,:);
        yt = yt+sar_res(1:n_nyq);
        % yt = yt+residual(1:n_nyq)
        yt = yt*(2/M/(M-1));

        fbin_sig=fsig/fnyq; w=hann_pv(n_nyq); fbin_L=3; fbin_H=n_nyq*fB/fnyq; N = n_nyq;
        [sqnr_test4, ptot_test4] = calcSNR(yt, fbin_sig, fbin_L, fbin_H, w, N);
        A = sqnr_test4
        sqnr_test4_array(1,end+1) = sqnr_test4;
    end

    figure(1+i_results); clf;
    hold on;
    scatter(Compen_Gain_Error_Array,sqnr_test3_array,[],'green');
    scatter(Compen_Gain_Error_Array,sqnr_test4_array,[],'magenta');
    peak_sqnr_test3 = max(sqnr_test3_array);
    peak_sqnr_test4 = max(sqnr_test4_array);
    i_peak_sqnr_test3 = find(sqnr_test3_array==peak_sqnr_test3);
    i_peak_sqnr_test4 = find(sqnr_test4_array==peak_sqnr_test4);
    value_compen_test3 = Compen_Gain_Error_Array(i_peak_sqnr_test3);
    value_compen_test4 = Compen_Gain_Error_Array(i_peak_sqnr_test4);
    text_handle= text(mean(Compen_Gain_Error_Array),min(sqnr_test3_array)+2, sprintf('Peak SQNR = %4.1fdB (x=%2d)',peak_sqnr_test3,value_compen_test3),'Color','green','FontSize',14,'HorizontalAlignment','center');
    text_handle= text(mean(Compen_Gain_Error_Array),min(sqnr_test3_array)+1, sprintf('Peak SQNR = %4.1fdB (x=%2d)',peak_sqnr_test4,value_compen_test4),'Color','magenta','FontSize',14,'HorizontalAlignment','center');
    legend('Division','Subtraction','Location','northwest');
    xlim([Compen_Gain_Error_Array(1),Compen_Gain_Error_Array(end)]);
    xlabel('x'); ylabel('SQNR [dB]');  
    title ({ ['Extended Counting',' (',num2str(SAR_B),'-Bit SAR, ',num2str(SAR_N_Comparison),' Quantization Cycles)']; ['V_P = ',num2str(sig_dB),' dBFS;  ','Fsig = ',num2str(fsig/1e3),' kHz;  ' ,'Fs = ',num2str(fs/1e6/OSR),' MHz; ','OSR = ', num2str(OSR)] });

    i_peak_sqnr_Devision_Array    = [i_peak_sqnr_Devision_Array    i_peak_sqnr_test3];
    i_peak_sqnr_Subtraction_Array = [i_peak_sqnr_Subtraction_Array i_peak_sqnr_test4];
end

%% C. DFT Analysis of (ExCount Dout + 10-B Quantized Analog Residue + Variable Coefficient:Compen_Gain_Error)
figure(11); clf;
hold on;
lengend_txt={};
for i_results=1:size(ptot_test_results,2)
    
    sar_res=D2_Array(i_results,:);
    if(Division_1__Subtraction_2==1)
        %%**** Division
        sar_res = (sar_res / (2^(SAR_B)+Compen_Gain_Error_Array(i_peak_sqnr_Devision_Array(i_results))));
    elseif(Division_1__Subtraction_2==2)
        %%**** Subtraction
        temp_sub = sar_res*Compen_Gain_Error_Array(i_peak_sqnr_Subtraction_Array(i_results)) / (2^(SAR_B));
        sar_res = ( (sar_res-temp_sub) / (2^(SAR_B)));
    else
        %%**** Default: No gain error compensation
        sar_res = (sar_res / (2^(SAR_B)));
    end
    yt = D1_Array(i_results,:);
    yt = yt+sar_res(1:n_nyq);
    yt = yt*(2/M/(M-1));

    [sqnr5_test, ptot5_test] = calcSNR(yt, fbin_sig, fbin_L, fbin_H, w, N);
    
    
    plot(linspace(0,fnyq/2,n_nyq/2+1), ptot5_test(1:n_nyq/2+1),'Color',Corner_color(i_results,:),'linewidth',2.0);
    test_message="";
    if(Division_1__Subtraction_2==1)
        %%**** Division
        text_handle =  text(5,-1-i_results*5, sprintf('SQNR = %5.1f dB (Coeff DEV=%d)',sqnr5_test,Compen_Gain_Error_Array(i_peak_sqnr_Devision_Array(i_results))),'Color',Corner_color(i_results,:),'FontSize',12);
    elseif(Division_1__Subtraction_2==2)
        %%**** Subtraction
%         text_handle =  text(5,-1-i_results*5, sprintf('SQNR = %5.1f dB (Coeff SUB=%d)',sqnr5_test,Compen_Gain_Error_Array(i_peak_sqnr_Subtraction_Array(i_results))),'Color',Corner_color(i_results,:),'FontSize',12);
        text_handle =  text(400e3,-1-i_results*5, sprintf('SQNR = %5.1f dB (Coeff SUB=%d)',sqnr5_test,Compen_Gain_Error_Array(i_peak_sqnr_Subtraction_Array(i_results))),'Color',Corner_color(i_results,:),'FontSize',12);
%         text_handle =  text(400e3,-1-i_results*5, sprintf('SQNR = %5.1f dB (Coeff SUB=%d)',sqnr5_test,Compen_Gain_Error_Array(i_peak_sqnr_Subtraction_Array(i_results))),'Color',Corner_color(i_results,:),'FontSize',12);
    else
        %%**** Default: No gain error compensation
        text_handle =  text(5,-1-i_results*5, sprintf('SQNR = %4.1f dB',sqnr5_test),'Color',Corner_color(i_results,:),'FontSize',12);
    end
%     text_handle= test_message;
    txt = [sprintf('%5s',num2str(Corner_PWRA(i_results))),'V ',Corner_Name{i_results},sprintf(' %3s',num2str(Corner_Temp(i_results))),char(0176),'C',];
    lengend_txt(end+1) = {txt};
end
xlabel('Frequency [Hz]'); ylabel('PSD [dB/Hz]');  axis([0 fB -130 0]); 
title ({['Extended Counting with 10-B Quantized Residue (',process_code,')'];['V_P= ',num2str(Vsig_dB),'dBFS;  ','Fsig=',num2str(fsig/1e3),'kHz;  ' ,'Fs=',num2str(fs/1e6/OSR),'MHz; ','OSR= ', num2str(OSR) ]},'FontSize',14);

% legend(lengend_txt);
legend(lengend_txt,'Location','northwest');
hold off;

%% D. DFT Analysis of (ExCount Dout + 10-B Quantized Analog Residue + Constant Coefficient:Compen_Gain_Error)
f12 = figure(12); clf;
set(gcf, 'color', [1 1 1]);
hold on;
lengend_txt={};
sqnr6_array=[];
for i_results=1:size(ptot_test_results,2)
    
    sar_res=D2_Array(i_results,:);
    if(Division_1__Subtraction_2==1)
        %%**** Division
        sar_res = (sar_res / (2^(SAR_B)+Compen_Gain_Error));
    elseif(Division_1__Subtraction_2==2)
        %%**** Subtraction
        temp_sub = sar_res*Compen_Gain_Error / (2^(SAR_B));
        sar_res = ( (sar_res-temp_sub) / (2^(SAR_B)));
    else
        %%**** Default: No gain error compensation
        sar_res = (sar_res / (2^(SAR_B)));
    end
    yt = D1_Array(i_results,:);
    yt = yt+sar_res(1:n_nyq);
    yt = yt*(2/M/(M-1));

    [sqnr6_test, ptot6_test] = calcSNR(yt, fbin_sig, fbin_L, fbin_H, w, N);
    sqnr6_array(end+1) = sqnr6_test;
    
    plot(linspace(0,fnyq/2/1e3,n_nyq/2+1), ptot6_test(1:n_nyq/2+1),'Color',Corner_color(i_results,:),'linewidth',2.0);
%     text_handle= text(5/1e3,-1-i_results*5, sprintf('SQNR = %4.1f dB',sqnr6_test),'Color',Corner_color(i_results,:),'FontSize',12);
%     text_handle= text(370e3/1e3,-1-i_results*5, sprintf('SQNR = %4.1f dB',sqnr6_test),'Color',Corner_color(i_results,:),'FontSize',12);
%         text_handle= text(480e3/1e3,-0-i_results*10.2, sprintf('SQNR = %4.1f dB',sqnr6_test),'Color',Corner_color(i_results,:),'FontSize',18);
%             text_handle= text(480e3/1e3,-0-i_results*10.2, sprintf('SQNR = %4.1f dB',sqnr6_test),'Color','k','FontSize',18);
%     text_handle= text(450e3/1e3,-1-i_results*5, sprintf('SQNR = %4.1f dB',sqnr6_test),'Color',Corner_color(i_results,:),'FontSize',12);
    txt = [sprintf('%5s',num2str(Corner_PWRA(i_results))),'V ',Corner_Name{i_results},sprintf(' %3s',num2str(Corner_Temp(i_results))),char(0176),'C',];
    lengend_txt(end+1) = {txt};
end
for i_results=1:size(ptot_test_results,2)
% 	text_handle= text(480e3/1e3,-0-i_results*10.2, sprintf('SQNR = %4.1f dB',sqnr6_array(i_results)),'Color',Corner_color(i_results,:),'FontSize',18);
	text_handle= text(480e3/1e3,-0-i_results*10.2, sprintf('SQNR = %4.1f dB',sqnr6_array(i_results)),'Color','k','FontSize',18);
end
xlabel('Frequency [kHz]','FontSize',32); ylabel('PSD [dB/Hz]','FontSize',32);  axis([0 fB/1e3 -130 0]); 
% title ({['Extended Counting with 10-B Quantized Residue (',process_code,')'];['V_P= ',num2str(Vsig_dB),'dBFS;  ','Fsig=',num2str(fsig/1e3),'kHz;  ' ,'Fs=',num2str(fs/1e6),'MHz; ','OSR= ', num2str(OSR),'; Coeff=',num2str(Compen_Gain_Error)]},'FontSize',14);
title ({['Extended Counting ADC (',process_code,')'];['V_P_P = ',num2str(Vsig_dB),' dBFS;  ','F_s_i_g = ',num2str(fsig/1e3),' kHz;  ' ,'F_B = ',num2str(fs/2/1e6/OSR),' MHz; ','OSR = ', num2str(OSR),'; x = ',num2str(Compen_Gain_Error)]},'FontSize',14);
set(gca, 'FontSize',13);
% scr_siz = get(0,'ScreenSize') ;
% scr_pos = f12.Position;
set(gcf, 'Position',  [488   342   580   420]);

% legend(lengend_txt);
legend(lengend_txt,'Location','northwest','FontSize',14);
legend('boxoff');
hold off;

%*********************************  END: Plotting Spectrum Analysis Results  ****************
%********************************************************************************************


