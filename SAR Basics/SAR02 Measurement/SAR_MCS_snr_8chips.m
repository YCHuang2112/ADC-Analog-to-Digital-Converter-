clear;

addpath ('D:\LocalLaptop\專題\研究所\adc\Matlab_Toolbox\schreier\delsig');
addpath ('D:\LocalLaptop\專題\研究所\adc\Matlab_Toolbox\[Malcovati] SDtoolbox\SDtoolbox');
% addpath ('D:\LocalLaptop\專題\研究所\adc\Matlab_Toolbox\MIT Design Tools\HSPICE_Toolbox\HspiceToolbox');
% addpath ('D:\LocalLaptop\專題\研究所\adc\Matlab_Toolbox\csvimport');

%% 8 Chips
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
MSB_CAL_ON = true; %(true,false)
% MSB_CAL_ON = false; %(true,false)

%% Set Input File & Signal Parameters
%% User's choice  <----------------------------------------------------------
filename_prefix = 'data\SAR02';
Chip_Dir = {'\chip1','\chip2','\chip3','\chip4','\chip5','\chip6','\chip7','\chip8'};
PWRA_USER=1.8; Vsig = 3.3; fsig=130e3; fsig_U=1e3; N_samples=1000; txt_sndr_xloc=0; txt_sndr_H='left'; txt_sfdr_xloc=0.18; txt_sfdr_H='left'; ConRate=1e6; ConRate_U=1e6;
% PWRA_USER=1.8; Vsig = 3.3; fsig=200e3; fsig_U=1e3; N_samples=1000; txt_sndr_xloc=0.01; txt_sndr_H='left'; txt_sfdr_xloc=0.21; txt_sfdr_H='left'; ConRate=1e6; ConRate_U=1e6;
% PWRA_USER=1.8; Vsig = 3.3; fsig=130e3; fsig_U=1e3; N_samples=1000; txt_sndr_xloc=0.01; txt_sndr_H='left'; txt_sfdr_xloc=0.4; txt_sfdr_H='left'; ConRate=2e6; ConRate_U=1e6;
% PWRA_USER=1.8; Vsig = 3.3; fsig=200e3; fsig_U=1e3; N_samples=1000;  txt_sndr_H='left'; txt_sfdr_H='left'; ConRate=2e6; ConRate_U=1e6; txt_sndr_xloc=0; txt_sfdr_xloc=0.38; %txt_sndr_xloc=0.2; txt_sfdr_xloc=0.5; 
if(fsig_U == 1e3)
    unit1 = 'k';
elseif(fsig_U == 1e6)
    unit1 = 'M';
end
if(ConRate_U == 1e3)
    unit2 = 'k';
elseif(ConRate_U == 1e6)
    unit2 = 'M';
end
filename_postfix= ['\',num2str(PWRA_USER),'V\',num2str(ConRate/ConRate_U),unit2,'Hz\sar02_mcs_meas_vin_',num2str(Vsig),'Vpp_',num2str(fsig/fsig_U),unit1,'Hz_ConRate',num2str(ConRate/ConRate_U),unit2,'Hz.csv'];

%% Set Power Supply Parameters
PWRA= PWRA_USER;
VFS = PWRA*2;

%% SAR-PART Parameters
N_bits = 10;  

%% File Lists
%%% SAR02　-----------------------------------------------------------------


%% set CLK Parameters
fs=ConRate; Ts=1/fs;  
fclk=fs*1; tclk=Ts;
fB = fs/2;

%% set Input-Signal Parameters
Vdb_sig = 20*log(Vsig/VFS);

%*********************************  END: Parameter settings  ********************************
%********************************************************************************************

%********************************************************************************************
%********************************* START: Post-Processing Loop  *****************************
sndr_test_results=[];sfdr_test_results=[];
enob_test_results=[];ptot_test_results=[];
for i_chip = 1:length(Chip_Dir)
    %********************************************************************************************
    %********************************** START: Signals Readin ***********************************
    %% Set File Name
    filename = [filename_prefix,Chip_Dir{i_chip},filename_postfix];
    
    %% Set Input Signal & Readin Signal
    for i=1:N_bits
        col_name = { char( "My Bus 1["+string(i-1)+"]" ) };
        readin = csvimport(filename,'columns',col_name);
        data(i,:) = readin';
    end
    
    %% Transform to Digital output
    digital_code = 0;
    for i=1:N_bits
        digital_code = digital_code + (data(i,:)-0.5)*2^(i-1);
    end
    if MSB_CAL_ON == true
            digital_code = digital_code - data(N_bits,:)*10/8;
%             for i_digital_code = 1:length(digital_code)
%                 if (digital_code(1,i_digital_code)*2) > 0
%                     digital_code(1,i_digital_code) = digital_code(1,i_digital_code)-10/16;
%                 elseif (digital_code(1,i_digital_code)*2) < 0
%                     digital_code(1,i_digital_code) = digital_code(1,i_digital_code)+10/16;
%                 end
%             end
    end
    digital_code = digital_code / 2^(N_bits);


    test = digital_code(end-N_samples+1-1200:end-1200);
%     test = digital_code(end-N_samples+1-2100:end-2100);
%     test = digital_code(end-N_samples+1-800:end-800);
    % test = digital_code(end-N_samples+1:end);


    %**********************************  END: Signals Readin  ***********************************
    %********************************************************************************************

    %********************************************************************************************
    %********************************* START: Spectrum Analysis *********************************
    %% DFT Analysis of SAR's DOUT
    fbin_sig_sp=fsig/fclk; w_sp=hann_pv(N_samples); fbin_L_sp=3; fbin_H_sp=N_samples*fB/fclk;
    yt=test;
    [sndr_test, ptot_test, a, b] = calcSNR(yt, fbin_sig_sp, fbin_L_sp, fbin_H_sp, w_sp, N_samples);

    shift = 7;
    fbin = fsig/fclk*N_samples+1;
    signal=(N_samples/sum(w_sp))*sinusx(yt(1:N_samples).*w_sp,fbin_sig_sp,N_samples);
    noise=yt(1:N_samples)-signal;	
    ssignal=(abs(fft((signal(1:N_samples).*w_sp)'))).^2;	% Signal PSD
    snoise=(abs(fft((noise(1:N_samples).*w_sp)'))).^2;		% Noise PSD
    pwsignal=sum(ssignal(fbin_L_sp:fbin_H_sp));	            % Signal power
    pwnoise=sum(snoise(fbin_L_sp:fbin-shift))+sum(snoise(fbin+shift:fbin_H_sp));		        % Noise power
    self_sndrdB=dbp(pwsignal/pwnoise);


    max_noise_db = max([ptot_test(fbin_L_sp:fbin-shift)' ptot_test(fbin+shift:fbin_H_sp)']);
    signal_db = ptot_test(floor(fbin));
    self_sfdrdB=signal_db-max_noise_db;
    
    sndr_test_results(end+1) = self_sndrdB;
    sfdr_test_results(end+1) = self_sfdrdB;
    enob_test_results(end+1) = ((self_sndrdB-1.76)/6.02);
    ptot_test_results(:,end+1) = ptot_test;
    %*********************************  END: Spectrum Analysis  *********************************
    %********************************************************************************************
end
%*********************************  END: Post-Processing Loop  ******************************
%********************************************************************************************

%********************************************************************************************
%********************************* START: Plotting Spectrum Analysis Results  ***************
%% Plot Analysis Results (DFT of SAR's DOUT)
f = 0:fs/N_samples:fs/2;
figure(1); clf;
set(gcf, 'color', [1 1 1]);
hold on;
lengend_txt={};
for i_chip=1:size(Chip_Dir,2)
    plot(f/ConRate_U,ptot_test_results(1:N_samples/2+1,i_chip),'Color',Corner_color(i_chip,:),'linewidth',2.0);
    plot(f/ConRate_U,10*log10(cumsum([ zeros(fbin_L_sp-1,1); 10.^(ptot_test_results(fbin_L_sp:fbin-shift,i_chip)/10); zeros(shift*2-1,1); 10.^(ptot_test_results(fbin+shift:fbin_H_sp+1,i_chip)/10); ])), 'Color',Corner_color(i_chip,:),'linewidth',2.0);
%     text_handle= text(floor(txt_xloc),-40, sprintf('ENOB　= %4.1fBits',enob_test_results(i_chip)),'Color','green','FontSize',18);
%     txt = [sprintf('%5s',num2str(Corner_PWRA(i_results))),'V ',Corner_Name{i_results},sprintf(' %3s',num2str(Corner_Temp(i_results))),char(0176),'C',];
    txt = ['chip',num2str(i_chip);];
    lengend_txt(end+1) = {txt};
end
for i_chip=1:size(Chip_Dir,2)
%     text_handle= text(txt_sndr_xloc,-2-i_chip*10, sprintf('SNDR　= %4.1fdB',sndr_test_results(i_chip)),'Color',Corner_color_txt(i_chip,:),'FontSize',12,'HorizontalAlignment',txt_sndr_H);
%     text_handle= text(txt_sfdr_xloc,-2-i_chip*10, sprintf('SFDR　= %4.1fdB',sfdr_test_results(i_chip)),'Color',Corner_color_txt(i_chip,:),'FontSize',12,'HorizontalAlignment',txt_sfdr_H);
    text_handle= text(txt_sndr_xloc,-2-i_chip*10, sprintf('%2d: SNDR　= %4.1fdB',i_chip,sndr_test_results(i_chip)),'Color','black','FontSize',12,'HorizontalAlignment',txt_sndr_H);
    text_handle= text(txt_sfdr_xloc,-2-i_chip*10, sprintf('SFDR　= %4.1fdB',sfdr_test_results(i_chip)),'Color','black','FontSize',12,'HorizontalAlignment',txt_sfdr_H);  
end
xlabel(['Frequency [',unit2,'Hz]'],'FontSize',20); ylabel('PSD [dB/Hz]','FontSize',20);  axis([0 fB/ConRate_U -120 0]); 
xlim([0 fB/ConRate_U]);

mean_sndr = mean(sndr_test_results);
sigma_sndr = sqrt(var(sndr_test_results));
title ({['10-bit SAR;  V_P_P = ',num2str(Vdb_sig),' dBFS;  '];['F_s_i_g = ',num2str(fsig/fsig_U),' ',unit1,'Hz;  ' ,'F_B = ',num2str(fB/ConRate_U),' ',unit2,'Hz;'];['MSB\_CAL\_ON = '+string(MSB_CAL_ON)]});
legend(lengend_txt,'FontSize',12);
legend('boxoff');
hold off;
set(gca, 'FontSize',14);
%*********************************  END: Plotting Spectrum Analysis Results  ****************
%********************************************************************************************

