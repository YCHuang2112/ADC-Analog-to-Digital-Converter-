function [Energy_per_switching, Energy_of_each_Dcode] = sar_MS_switching(nth_bit,B_SAR,i_index,DACP_pre,DACN_pre,DACP_cur,DACN_cur,Energy_per_switching,Q_CA_CB_RedistributionTable_MS,DigitalCode,E_cur,Energy_of_each_Dcode)
    if(nth_bit == B_SAR)
%         disp(DigitalCode);
%         disp(E_cur);
        Energy_of_each_Dcode(2*DigitalCode+1) = E_cur;
        Energy_of_each_Dcode(2*DigitalCode+1+1) = E_cur;
        return;
    end
%     fprintf(' nth_bit=%d i_index=%d \n',nth_bit, i_index); 
%     disp(DACP_cur);
%     disp(DACN_cur);
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
    DACP_next = DACP_cur;    DACN_next = DACN_cur;
    DACN_next(nth_bit)=0;   
    
    %%% Switching Energy on CDACp
    if(nth_bit < B_SAR)
        dVn1 = -1*Q_CA_CB_RedistributionTable_MS(nth_bit,:);
        En = dVn1*[DACN_next'];
        Etot = En;
        Energy_per_switching(i_index*2)=Etot;
    end
    
    [Energy_per_switching,Energy_of_each_Dcode] = sar_MS_switching(nth_bit+1,B_SAR,i_index*2,DACP_cur,DACN_cur,DACP_next,DACN_next,Energy_per_switching,Q_CA_CB_RedistributionTable_MS,DigitalCode*2,E_cur+Etot,Energy_of_each_Dcode);

    
    %%% ANS: Yes
    DACP_next = DACP_cur;    DACN_next = DACN_cur;
    DACP_next(nth_bit)=0;
    
    %%% Switching Energy on CDACp
    if(nth_bit < B_SAR)
        dVp1 = -1*Q_CA_CB_RedistributionTable_MS(nth_bit,:);
        Ep = dVp1*[DACP_next'];
        Etot = Ep;
        Energy_per_switching(i_index*2+1)=Etot;
    end
    
    [Energy_per_switching,Energy_of_each_Dcode] = sar_MS_switching(nth_bit+1,B_SAR,i_index*2+1,DACP_cur,DACN_cur,DACP_next,DACN_next,Energy_per_switching,Q_CA_CB_RedistributionTable_MS,DigitalCode*2+1,E_cur+Etot,Energy_of_each_Dcode);

end
