function [Energy_per_switching, Energy_of_each_Dcode] = sar_MCS_switching(nth_bit,B_SAR,i_index,DACP_Vcm_pre,DACN_Vcm_pre,DACP_Vcm_cur,DACN_Vcm_cur,DACP_Vref_pre,DACN_Vref_pre,DACP_Vref_cur,DACN_Vref_cur,Energy_per_switching,Q_CA_CB_RedistributionTable_MCS,DigitalCode,E_cur,Energy_of_each_Dcode)
    if(nth_bit == B_SAR)
        disp(DigitalCode);
        disp(E_cur);
        Energy_of_each_Dcode(2*DigitalCode+1) = E_cur;
        Energy_of_each_Dcode(2*DigitalCode+1+1) = E_cur;
        return;
    end
%     fprintf(' nth_bit=%d i_index=%d \n',nth_bit, i_index); 
%     disp(DACP_Vref_cur);
%     disp(DACN_Vref_cur);
%     for i=1:B_SAR+1 
%         fprintf('%d ',DACP_pre(i));
%     end
%     fprintf("\n");
%     for i=1:B_SAR+1
%         fprintf('%d ',DACN_pre(i));
%     end
%     fprintf("\n\n");
    
    %Q: V(DACP) > V(DACN)?
    %%% ANS: No
    DACP_Vref_next = DACP_Vref_cur;    DACN_Vref_next = DACN_Vref_cur;
    DACP_Vref_next(nth_bit)=1;    
    DACN_Vref_next(nth_bit)=0;    
    
    DACP_Vcm_next = DACP_Vcm_cur;    DACN_Vcm_next = DACN_Vcm_cur;
    DACP_Vcm_next(nth_bit)=0;    
    DACN_Vcm_next(nth_bit)=0;    
    
    %%% Switching Energy on CDACp
    if(nth_bit < B_SAR)
        dVp1 = 0.5*Q_CA_CB_RedistributionTable_MCS(nth_bit,:);
        Ep = dVp1*[DACP_Vref_next'];
        dVn1 = -0.5*Q_CA_CB_RedistributionTable_MCS(nth_bit,:);
        En = dVn1*[DACN_Vref_next'];
        Etot = Ep + En;
        Energy_per_switching(i_index*2)=Etot;
    end
    
    [Energy_per_switching,Energy_of_each_Dcode] = sar_MCS_switching(nth_bit+1,B_SAR,i_index*2,DACP_Vcm_cur,DACN_Vcm_cur,DACP_Vcm_next,DACN_Vcm_next,DACP_Vref_cur,DACN_Vref_cur,DACP_Vref_next,DACN_Vref_next,Energy_per_switching,Q_CA_CB_RedistributionTable_MCS,DigitalCode*2,E_cur+Etot,Energy_of_each_Dcode);

    
    %%% ANS: Yes
    DACP_Vref_next = DACP_Vref_cur;    DACN_next = DACN_Vref_cur;
    DACP_Vref_next(nth_bit)=0;
    DACN_Vref_next(nth_bit)=1;   
    
    DACP_Vcm_next = DACP_Vcm_cur;    DACN_Vcm_next = DACN_Vcm_cur;
    DACP_Vcm_next(nth_bit)=0;    
    DACN_Vcm_next(nth_bit)=0;    
    
    %%% Switching Energy on CDACp
    if(nth_bit < B_SAR)
        dVp1 = -0.5*Q_CA_CB_RedistributionTable_MCS(nth_bit,:);
        Ep = dVp1*[DACP_Vref_next'];
        dVn1 = 0.5*Q_CA_CB_RedistributionTable_MCS(nth_bit,:);
        En = dVn1*[DACN_Vref_next'];
        Etot = Ep + En;
        Energy_per_switching(i_index*2+1)=Etot;
    end
    
    [Energy_per_switching,Energy_of_each_Dcode] = sar_MCS_switching(nth_bit+1,B_SAR,i_index*2+1,DACP_Vcm_cur,DACN_Vcm_cur,DACP_Vcm_next,DACN_Vcm_next,DACP_Vref_cur,DACN_Vref_cur,DACP_Vref_next,DACN_Vref_next,Energy_per_switching,Q_CA_CB_RedistributionTable_MCS,DigitalCode*2+1,E_cur+Etot,Energy_of_each_Dcode);

end
