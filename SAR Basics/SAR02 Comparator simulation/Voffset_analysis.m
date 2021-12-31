clear all;
 
%%%%%%%%%%% parameters %%%%%%%%%%%
vdd = 1.8; vss=0; vcm=0.9; num_mc=50; point_in_lsb=100;
vp=0.9; vfs=vp*2; n_bits=10; vlsb = vfs/(2^n_bits);
%%%%%%%%%%% parameters %%%%%%%%%%%
%%%% Comparator
N_LSB = 8;
filename = 'Voffset_Fs100M_8LSB_p100.m';

run(filename);

test = v__outp__(:,1:num_mc);
test_i = find(test > vcm);
test(test_i) = 1;
test_i = find(test < vcm);
test(test_i) = 0;
dis = sum(test,2);
dis = dis(2:end);
% accu = dis;
dis = dis';
len_dis = length(dis);
dis = [dis(1:(len_dis+1)/2) num_mc-dis((len_dis+1)/2+1:end)];
% plot(dis);

tmp_array = [];
for i=1:len_dis
    for j=1:dis(i)
        tmp_array(end+1) = i;
    end
end



vov_std_in_lsb   = std(tmp_array/point_in_lsb);
vov_std_in_lsb_v = vov_std_in_lsb * vlsb;
vov_3std_in_lsb   = vov_std_in_lsb *3
vov_3std_in_lsb_v = vov_std_in_lsb_v *3

figure(1);
set(gcf, 'color', [1 1 1]);
plot([-N_LSB/2*point_in_lsb:N_LSB/2*point_in_lsb]/point_in_lsb,dis);
title("Comparator Offset",'FontSize',32);
text_handle= text(0,28, sprintf('3 \\sigma = %4.1f LSB',vov_3std_in_lsb),'Color','green','FontSize',32,'HorizontalAlignment','center');
xlabel("(LSB)",'FontSize',32);
ylabel("Count(s)",'FontSize',32);
set(gca, 'FontSize',20);