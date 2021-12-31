clear;
Lth_order = 5;
M_cycles = 128;
i_start_M_cycles=1;
i_end_M_cycles=128;
sum = [];
square_of_sum = [];
sum_of_square = [];
Gth_L = [];


% i_order == 1
sum = ones(1,M_cycles);

% i_order > 1
for i_order = 2:(Lth_order+1)
    for i_cycle = 1:M_cycles
        if(i_cycle < i_order)
            sum(i_order,i_cycle)=0;
        else
            sum(i_order,i_cycle)=sum(i_order,i_cycle-1)+sum(i_order-1,i_cycle-1);
        end
    end
end

square_of_sum = sum.*sum;

sum_of_square = zeros(Lth_order,1);
for i_order = 1:Lth_order
    for i_cycle = 2:M_cycles
        if(i_cycle < i_order)
            sum_of_square(i_order,i_cycle)=0;
        else
            sum_of_square(i_order,i_cycle)=sum_of_square(i_order,i_cycle-1)+square_of_sum(i_order,i_cycle-1);
        end
    end
end

square_of_sum(square_of_sum==0)=-1;
M_by_sum_of_square = [1:M_cycles].*sum_of_square;
Gth_L = M_by_sum_of_square./square_of_sum(2:end,1:end);

figure(1); %%*** With delay Integrators
set(gcf, 'color', [1 1 1]);
Gth_L = Gth_L(:,i_start_M_cycles:i_end_M_cycles);
plot(Gth_L');



sum = [];
square_of_sum = [];
sum_of_square = [];
Gth_L = [];


% i_order == 1
sum = ones(Lth_order+1,M_cycles);

% i_order > 1
for i_order = 2:(Lth_order+1)
    for i_cycle = 2:M_cycles
        sum(i_order,i_cycle)=sum(i_order,i_cycle-1)+sum(i_order-1,i_cycle);
    end
end

square_of_sum = sum.*sum;

sum_of_square = ones(Lth_order,1);
for i_order = 1:Lth_order
    for i_cycle = 2:M_cycles
        sum_of_square(i_order,i_cycle)=sum_of_square(i_order,i_cycle-1)+square_of_sum(i_order,i_cycle);
    end
end

% square_of_sum(square_of_sum==0)=-1;
M_by_sum_of_square = [1:M_cycles].*sum_of_square;
Gth_L = M_by_sum_of_square./square_of_sum(2:end,1:end);

figure(2); %%*** With non-delay integrators
set(gcf, 'color', [1 1 1]);
Gth_L = Gth_L(:,i_start_M_cycles:i_end_M_cycles);
for i_order = Lth_order:-1:1
    data(:,Lth_order-i_order+1) = Gth_L(i_order,:);
end
semilogx(data,'LineWidth',2);
xlim([1 128]);
ylim([0.5 3]);
x_axis_tick = 2.^(0:7);
xticks(x_axis_tick);
set(gca, 'FontSize', 16);
xlabel('M (Cycles)','FontSize',18);
ylabel('G_t_h_,_L(M)','FontSize',18);
legend('L=5','L=4','L=3','L=2','L=1',"location",'northwest','FontSize',16);
% title('thermal noise penalty');