clc,clear

%==========速率方程不考虑  噪声 相位  =============================================
%==========对应方程 rate_equ_nophase =======


%% 参数

q = 1.6e-19;    %C,电子电量
c = 3e10;       %光速度，单位（cm/s）
h = 6.62e-34; 	%普朗克常数

V = 4e-12;      %有源区体积 cm^3
conf = 0.032;     %光限制因子
ng = 4.2;         %群折射率
vg = c/ng;        %群速度 cm/s

g0 = 1800;        %线性增益系数   cm^-1
Ntr = 1.8e18;     %透明载流子浓度 cm^-3
Ns = -0.4e18;     %增益参数       cm^-3
eps = 1.5e-17;    %增益压缩因子    cm^3

A = 0;            %非辐射复合系数 
B = 0.8e-10;      %辐射复合系数  cm^3/s
C = 3.5e-30;      %俄歇复合系数  cm^6/s

alpha = 5;        %线宽增益因子
eta = 0.8;        %电流注入效率
eta_0 = 0.45;     %光收集效率
beta = 0.895e-4;  %自发发射因子

lambda = 0.98e-4;	%波长，cm
mu = c/lambda;    %角频率，Hz
tp = 2.77e-12;    %光子寿命,s
a_i = 5;        %内损耗，cm-1
a_m = 45.6;      %镜面损耗，cm-1

%% Part1 直流偏置下，载流子,光子数 随电流的变化,可得出 阈值电流，阈值载流子，阈值增益


t1 = 0; 
t2 = 5e-9; 
dt = 5e-14; 
nsteps = (t2-t1)/(dt) +1;
tspan = linspace(t1,t2,nsteps); %span of time values

Im = 1e-3;            %调制幅度
f = 2e9;              %调制频率，Hz
on = 0;               %调制=0
I=0:0.1e-3:2e-3;      %电流0-20mA
% y0 = [1.8e18;0];
y0 = [10;10];

N_i = zeros(1,length(I));
delta_P_w = zeros(1,length(I));
delta_Pw = zeros(1,length(delta_P_w));
gain_i = zeros(1,length(I));

for i=1:1:length(I)
    I0 = I(i);
    for j=1:3
    [t,y] = ode23(@(t,y) rate_equ_nophase(t,y,I0,Im,f,on), tspan, y0);  %计算结果有虚数，换成odde23未能解决
    y0=[y(length(tspan),1);y(length(tspan),2)];
    end
    L=length(y); 
    N_i(i) = y(L,1);          %不同电流下稳态载流子
    delta_P_w(i) = y(L,2);    %不同点留下稳态光子数
end

delta_Pw = eta_0*h*mu.*delta_P_w*(V/conf)*vg*a_m;         	%功率
gain_i = g0./(1+eps.*delta_P_w).*log((N_i+Ns)/(Ntr+Ns));      %增益

figure                        %载流子,光子数 随电流的变化

yyaxis left
plot(1000*I,N_i/(1e18),'-o')
text(0.6,1.5,'g')
hold on 
plot(1000*I,1000*delta_Pw*10,'-d','color','red')
text(0.5,3.2,'N')
xlabel('Current(mA)','Fontsize',10,'color','black');
ylabel('N(\times 10^{18} cm^{-3}),  P(\times 0.1mW）','Fontsize',10,'color','black');

yyaxis right
plot(1000*I,gain_i,'-*')
text(1.5,800,'P')
ylim([0 2000])
ylabel('Gain(cm^{-1})','Fontsize',10,'color','black');
hold off
axis square
% 
% 
% 
% figure
% plot(N_i,gain_i,'-o'); 	%增益随载流子的变化
% xlabel('Carrier density (cm^{-3})','Fontsize',10,'color','black');
% ylabel('gain(cm^{-1})','Fontsize',10,'color','black');
% grid on
% axis square


%% Part2 直流偏置或者调制状态下 载流子与光子随时间变化关系

t1 = 10e-9;     	%Start time  %s微分方程求解时间区间,从td+t1->td+t2
t2 = 100e-9;        %Stop time
dt = 5e-12;         %步长
nsteps = (t2-t1)/(dt) +1;
tspan = linspace(t1,t2,nsteps);  


f = 2e9;            %调制频率，Hz
I0 = 3e-3;          %偏置电流
Im = 2e-3;          %调制幅度
on = 0             %调制=0



y0 = [1.8e18;1e13];                                                 
[t,y] = ode23(@(t,y) rate_equ(t,y,I0,Im,f,on),tspan, y0);   
I = I0+Im*sin(2*pi*f*t)*on;
    for j=1:3
    [t,y] = ode23(@(t,y) rate_equ_nophase(t,y,I0,Im,f), tspan, y0);  %计算结果有虚数，换成odde23未能解决
    y0=[y(length(tspan),1);y(length(tspan),2)];
    end


figure

yyaxis left
plot(t*1e9,y(:,1))             %ns-N/1e18
% hold on 
% plot(t*1e9,I*1e3,'color','red')     %ns-mA
xlim([10 15])
hold off
% ylim([0 6])

xlabel('t(ns)','Fontsize',10,'color','black');
ylabel('Carrier density','Fontsize',10,'color','black');

yyaxis right
plot(t*1e9,y(:,2))
ylabel('Photon density','Fontsize',10,'color','black');

axis square

%% Part3 解析解求小信号调制参数

%直流偏置,特定电流下获取稳态光子数,小信号参数================================


% f = 2e9;            %调制频率，Hz
% I0 = 3e-3;          %偏置电流
% Im = I0-1e-3;       %调制幅度
% 
% 
% on = 0;                       %调制=0
% 
% t1 = 1e-9; % Start time       %微分方程求解时间区间,从td+t1->td+t2
% t2 = 5e-9; % Stop time
% dt = 5e-14;                   %步长
% nsteps = (t2-t1)/(dt) +1;
% tspan = linspace(t1,t2,nsteps); %span of time values
% y0 = [1e18;0];
% 
% % for j=1:3
% % [t,y] = ode45(@(t,y) rate_equ(t,y,I0,Im,f,on), tspan, y0);  %计算结果有虚数，换成odde23未能解决
% % y0=[y(length(tspan),1);y(length(tspan),2)];
% % end
% [t,y] = ode45(@(t,y) rate_equ(t,y,I0,Im,f,on), tspan, y0);
% I = I0+Im*sin(2*pi*f*t)*on;
% 
% L=length(y);
% N = y(L,1)                    %稳态电子数
% P = y(L,2)                    %稳态光子数
% figure
% plot(t,y(:,2))                %光子数稳定过程
% 
% %利用N,P值，求小信号参数===================================================
% 
% gain = rate_equ_gain(N,P);
% a0 = 5.3386e-16;              %cm^2
% ap = eps*gain/(1+eps*P);
% a = a0/(1+eps*P);
% 
% rnn = A+2*B*N+3*C*N^2+P*vg*a;
% rnp = 1/conf/tp-vg*ap*P;
% rpn = conf*vg*a*P;
% rpp = conf*vg*ap*P;
% 
% wr = (rnp*rpn+rnn*rpp)^0.5;
% fr = wr/(2*pi)  %Hz
% r = rnn+rpp     %Hz
% 
% f = 0:0.01e9:10e9;            %0-10GHz调制响应
% H = zeros(1,length(f));
% 
% H = wr^2./(wr^2-(2*pi*f).^2+1i*2*pi*f*r);
% Hw = 10*log10(abs(H).^2);
% 
% figure
% plot(f/1e9,Hw)                %小信号响应曲线 横坐标GHz
% xlabel('Frequency(GHz)','Fontsize',10,'color','black');
% ylabel('Respones','Fontsize',10,'color','black');
% grid on
% axis square


%% Part4 速率方程解小信号响应

% % % I0 = 3e-3;%偏置电流
% % % Im = 0.5e-3;%调制幅度
% I = I0+Im*sin(2*pi*f*t)*on;

% % % on = 1;                       %调制=1
% % % 
% % % 
% % % mod_f=0.2e9:0.1e9:10e9;      %调制频率
% % % 
% % % t1 = 1e-9; 
% % % t2 = 50e-9; 
% % % dt = 5e-13; 
% % % nsteps = (t2-t1)/(dt) +1;
% % % tspan = linspace(t1,t2,nsteps); 
% y0 = [1.8e18;0];
% % % y0 = [10;10];
% % % 
% % % 
% % % 
% % % 
% % % delta_P_w = zeros(1,length(mod_f));
% % % delta_Pw = zeros(1,length(delta_P_w));


% for i=1:1:length(mod_f)
% % %     mod_f(50)
% % %     f = mod_f(50);
% % %  
% % %             [t,y] = ode45(@(t,y) rate_equ(t,y,I0,Im,f,0), tspan, y0);  
% % %             y0=[y(length(tspan),1);y(length(tspan),2)];
% % %             [t,y] = ode45(@(t,y) rate_equ(t,y,I0,Im,f,on), tspan, y0);

%     ph = y(end-8000:end,2);
%     f_max = find(localmaxmin(ph,'max')); %极大值的位置
%     jmax = ph(f_max);
%     f_min = find(localmaxmin(ph,'min'));
%     jmin = ph(f_min);
% 
%     delta_P_w(i) = 0.5*(  mean(jmax) - mean(jmin) );       %不同调制频率下，光子数变化的幅度
% end


% Hw = delta_P_w/Im/(eta*conf*tp/q/V);                    %不同调制频率下，调制响应,书中公式（5.56）
% R = 20*log10(abs(Hw));                                  %不同调制频率下，调制响应
% 
% figure(1)               %调制电流曲线，对应的光子数变化曲线，此时调制频率为 max(mode_f)

% % % I = I0+Im*sin(2*pi*f*t)*on;
% % % 
% % % yyaxis left
% % % plot(t,I)
% % % % xlim([0 1e-9])
% % % xlabel('t','Fontsize',10,'color','black');
% % % ylabel('Current','Fontsize',10,'color','black');
% % % 
% % % yyaxis right
% % % plot(t,y(:,2))
% % % ylabel('Photon','Fontsize',10,'color','black');
% % % 
% % % axis square


% figure(2)               %调制响应曲线
% plot(mod_f/1e9,R,'-o');     
% xlabel('Frequency(GHz)','Fontsize',10,'color','black');
% ylabel('Response','Fontsize',10,'color','black');
% grid on
% axis square
% plot(mod_f/1e9,delta_P_w,'-o'); 

%%


















