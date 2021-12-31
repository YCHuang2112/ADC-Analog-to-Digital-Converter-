clear;

addpath ('D:\LocalLaptop\專題\研究所\adc\Matlab_Toolbox\csvimport');

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
N_bits = 10;
MSB_CAL_ON = true; %(true,fatlse)
% MSB_CAL_ON = false; %(true,false)

%% Set Input File & Signal Parameters
%% User's choice  <----------------------------------------------------------
filename_prefix = 'data\SAR02\INLDNL';
Chip_Dir = {'\chip1','\chip2','\chip3','\chip4','\chip5','\chip6','\chip7','\chip8'};
PWRA_USER=1.8; Vsig = 3.5; fsig=2.04e3; fsig_U=1e3; N_samples=600; txt_sndr_xloc=0.5; txt_sndr_H='left'; txt_sfdr_xloc=1.5; txt_sfdr_H='left'; ConRate=10e3; ConRate_U=1e3;
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
filename_postfix= ['\vin_',num2str(Vsig),'Vpp_',num2str(fsig/fsig_U),unit1,'Hz_ConRate',num2str(ConRate/ConRate_U),unit2,'Hz_Np4M_Dout_Hist.csv'];

%% File Lists
%%% SAR02　-----------------------------------------------------------------
%%% **  CLK--SILICON LABS Si5338-EVB REV 1.0, SIGNAL--Audio Precision SYS-2722

%*********************************  END: Parameter settings  ********************************
%********************************************************************************************

%********************************************************************************************
%********************************* START: 8-Chip Loop  **************************************
sndr_test_results=[];sfdr_test_results=[];
enob_test_results=[];ptot_test_results=[];
Gain_Error_results=[];
for i_chip = 1:length(Chip_Dir)
    %********************************************************************************************
    %********************************** START: Signals Readin ***********************************
   %% Set File Name
    filename = [filename_prefix,Chip_Dir{i_chip},filename_postfix];
    %% Set Input Signal & Readin Signal

    col_name='Dout_Hist';
    hv = csvimport(filename,'columns',col_name);
    hv = hv';
    %**********************************  END: Signals Readin  ***********************************
    %********************************************************************************************

    %********************************************************************************************
    %********************************** START: Compute INL/DNL **********************************
    % hv = h.Values;
    % hv(1:1) = 0;
    % hv(1024:1024) = 0;
    cum = cumsum(hv);
    % figure(1);
    % plot(hv);

    % 
    % figure(2);
    % plot(cum);

    ch = sin((cum/max(cum)-0.5)*pi);
    % figure(3);
    % plot(ch);


    bl=1; bh=1022;
    charac_line = (ch+1)/2*(bh-bl+1);
    % figure(4);
    % hold on;
    % plot(charac_line);
    % % plot([1:1024],[1:1024]-47);
    % hold off;

    % dnl = [1 charac_line(2:1024)-charac_line(1:1023)]-1;
    % figure(5);
    % plot(dnl);
    % ylabel("(LSB)");

%     inl = cumsum(dnl(bl:bh));
    % inl = cumsum(dnl);
    m = [1:1024]'; n = ones(1024,1);
    if MSB_CAL_ON == true
        m_cal = [[1:511]+10/16,512,[513:1024]-10/16]; m=m_cal';
    end
    A = [m n];
    y = charac_line';
    x = A\y;
    
    % figure(6);
    % hold on;
    % plot(1:1024,y,'r',1:1024,A*x,'g');
    % hold off;

    inl2 =  y-A*x;
    inl2 = inl2';
    dnl2 = [0 inl2(2:1024)-inl2(1:1023)];
    INL2_ARR(i_chip,:) = inl2;
    DNL2_ARR(i_chip,:) = dnl2;
    
    Gain_Error_results(end+1)=1/x(1);
    %**********************************  END: Compute INL/DNL  **********************************
    %********************************************************************************************
end
%**********************************  END: 8-Chips Loop  *************************************
%********************************************************************************************
    
%********************************************************************************************
%********************************** START: Plotting INL/DNL Analysis Results ****************
figure(7);
set(gcf, 'color', [1 1 1]);
hold on;
lengend_txt={};
for i_chip = 1:length(Chip_Dir)
    % plot(inl2(1:1023));
    plot(INL2_ARR(i_chip,bl:bh));
    txt = ['chip',num2str(i_chip),' (gain error = ',num2str(Gain_Error_results(i_chip)),')';];
    lengend_txt(end+1) = {txt};
end
hold off;
% title(['INL; Gain Error = ',num2str(x(1))]);
title(['INL; MSB\_CAL\_ON = '+string(MSB_CAL_ON)],'FontSize',32);
% legend(lengend_txt,'Location','southwest','FontSize',16,'NumColumns',1','color','none');
legend(lengend_txt,'Location','southwest','FontSize',12,'NumColumns',2','color','none');
legend('boxoff');
xlim([0 1023]);
ylim([-2.5 1.5]);
xlabel('Digital Code','FontSize',20); ylabel('LSB','FontSize',20);
set(gca, 'FontSize',20);


figure(8);
set(gcf, 'color', [1 1 1]);
hold on;
lengend_txt={};
for i_chip = 1:length(Chip_Dir)
    % plot(dnl2(1:1023));
    plot(DNL2_ARR(i_chip,bl:bh));
    txt = ['chip',num2str(i_chip),' (gain error = ',num2str(Gain_Error_results(i_chip)),')';];
    lengend_txt(end+1) = {txt};
end
hold off;
% title(['DNL; Gain Error = ',num2str(x(1))]);
title(['DNL; MSB\_CAL\_ON = '+string(MSB_CAL_ON)],'FontSize',32);
legend(lengend_txt,'Location','southwest','FontSize',12,'NumColumns',2);
legend('boxoff');
xlim([0 1023]);
ylim([-1 0.6]);
xlabel('Digital Code','FontSize',20); ylabel('LSB','FontSize',20);
set(gca, 'FontSize',20);
%**********************************  END: Plotting INL/DNL Analysis Results  ****************
%********************************************************************************************


