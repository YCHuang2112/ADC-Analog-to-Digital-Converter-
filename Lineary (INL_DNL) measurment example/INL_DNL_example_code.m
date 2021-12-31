B = 8;
n_bits = 8;
range = 2^(n_bits-1) - 1;
vrefp = range           
vrefn = -range


% vdelta = (vrefp-vrefn)/(2^B-1);
% range = 2^(n_bits-1) - 1;
th = vrefn: vrefp;
 th(20) = th(20) + 0.7;

fs = 10e6;
fx = 49.4e6+314;
% fx = 494e3;
C = round(100 * 2^B / (fs/fx));

t = 0:1/fs:C/fs;
%sine input
x = (range+1) * sin(2*pi*fx.*t);
% ramp input
% x = (range+1) * (t/C*fs*2-1);

n_input = length(x);
y = zeros(1,n_input)
for i = 1:n_input
    for j=1:(2^n_bits-1)
%         vlev  =  vrefn + vdelta * (j-1);
%         vfloor = vlev - vdelta/2;
%         vcell  = vlev + vdelta/2;
%         if( (vfloor <= x(i))&&(x(i)<vcell) ) 
%             y(i) = (j-1);
%         end
          if(x(i) > th(j))
              y(i) = j;
          end
    end
end

y = y - 2^(n_bits-1);

figure(1);
hist(y,min(y):max(y));

% dnl_inl_sin(y);
minbin=min(y);
maxbin=max(y);

h = hist(y, minbin:maxbin);
ch = cumsum(h);
figure(6);
plot(pi*ch/sum(h)-pi/2);
figure(7);
plot(sin(pi*ch/sum(h)-pi/2));
T = -cos(pi*ch/sum(h));
% figure(2);
% plot(ch);
hlin = T(2:end) - T(1:end-1);
trunc = 2;
hlin_trunc = hlin(1+trunc:end-trunc);

lsb = sum(hlin_trunc)  / (length(hlin_trunc));
dnl  = [0 hlin_trunc/lsb-1];
misscodes = length(find(dnl<-0.9));

inl = cumsum(dnl);
% 
figure(8);
plot(dnl,'o');
ylim([-2 2]);
title("dnl");
figure(9);
plot(inl,'o');
ylim([-2 2]); 
title("inl");



figure(10);
x = -pi/2:pi/100:pi/2;
% sine input
% plot(sin(x)*(range+1));
% ramp input
plot(x/pi*2*(range+1));
line([0,100],[128,128],'Color',[17 17 17]/255,'LineStyle','--');
line([0,100],[-128,-128],'Color',[17 17 17]/255,'LineStyle','--');
for i = -127:127
    line([105,110],[i,i],'Color',[255 0 0 80]/255);
end
line([105,110],[128,128],'Color',[17 17 17]/255,'LineStyle','--');
line([105,110],[-128,-128],'Color',[17 17 17]/255,'LineStyle','--');