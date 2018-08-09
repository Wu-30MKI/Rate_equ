%% 速率方程
%====== 不包含噪声，不包括外光注入,包括相位===对应方程 rate_equ_nonoise==============
%====== 仿真 无外光注入下的 静态特性 和 动态特性 ==========================
%====== 通过改变 自发辐射因子beta，得到主激光器的参数 Sm Nm Phim ======================
% clc 
% clear
% clf
% close all
%% 激光器参数
q = 1.6e-19;    %C,电子电量
c = 3e10;       %光速度，单位（cm/s）
h = 6.62e-34;  	%普朗克常数
V = 4e-12;      %有源区体积 cm^3
conf=0.032;     %光限制因子
ng=4.2;         %群折射率
vg=c/ng;      	%群速度 cm/s

g0=1800;        %线性增益系数 cm^-1
Ntr=1.8e18;     %透明载流子浓度 cm^-3
Ns=-0.4e18;     %增益参数 cm^-3
eps=1.5e-17;    %增益压缩因子 cm^3

A=0;            %非辐射复合系数 
B=0.8e-10;      %辐射复合系数  cm^3/s
C=3.5e-30;      %俄歇复合系数  cm^6/s
alpha=5;        %线宽增益因子
eta=0.8;        %电流注入效率
beta=0.895e-6;      %自发发射因子

lambda=0.98e-4;	%波长，cm
mu=c/lambda;    %角频率，Hz
tp=2.77e-12;    %光子寿命,s   tau_p

%% 直流偏置下，载流子,光子数 随电流的变化,可得出 阈值电流，阈值载流子，阈值增益

t1 = 0; 
t2 = 5e-9; 
dt = 5e-14;
nsteps = (t2-t1)/(dt) +1;
tspan = linspace(t1,t2,nsteps); 
Bsim = 1/dt;

Im = 1e-3;            %调制幅度
f = 2e9;              %调制频率，Hz
on = 0;               %调制=0
I=0:0.1e-3:1e-3;      %电流0-10mA
% y0 = [1.8e18;1e13;1]; 
y0 = [Ntr;1e15;1];

N_i = zeros(1,length(I));           %载流子数，
P_i = zeros(1,length(I));           %光子数，
Power = zeros(1,length(P_i));       %功率，W
gain_i = zeros(1,length(I));

for i=1:1:length(I)
    I0 = I(i);
    [t,y] = ode23(@(t,y) rate_equ_nonoise(t,y,I0,Im,f,on,Bsim),tspan, y0);
%     for j=1:3                                            
%     [t,y] = ode23(@(t,y) rate_equ_switch(t,y,I0,Im,f,on,Bsim),tspan, y0);
%     y0=[y(length(tspan),1);y(length(tspan),2);y(length(tspan))];
%     end
    L=length(y); 
    N_i(i) = y(L,1);          %不同电流下稳态载流子
    P_i(i) = y(L,2);    %不同电流下稳态光子数
end

Power = eta_0*h*mu.*P_i*(V/conf)*vg*a_m;               %功率
% % gain_i = g0./(1+eps.*delta_P_w).*log((N_i+Ns)/(Ntr+Ns));        %增益


figure            %直流偏置下，载流子,光子数随电流的变化

yyaxis left
plot(1000*I,N_i/(1e18),'-o')
text(0.6,1.5,'g')
hold on 
plot(1000*I,1000*Power*10,'-d','color','red')
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



% figure(2)                 %增益随载流子的变化
% plot(N_i,gain_i,'-o'); 	
% xlabel('Carrier density (m^(-3))','Fontsize',10,'color','black');
% ylabel('gain(cm^{-1})','Fontsize',10,'color','black');
% grid on
% axis square



%% 直流偏置或者调制状态下 载流子与光子随时间变化关系

t1 = 10e-9;     	%Start time  %s微分方程求解时间区间,从td+t1->td+t2
t2 = 50e-9;        %Stop time
dt = 1e-12;         %步长
nsteps = (t2-t1)/(dt) +1;
tspan = linspace(t1,t2,nsteps);    

Im = 1e-3;       %调制幅度
f = 5e9;            %调制频率，Hz
I0 = 3e-3;          %偏置电流
on = 0;             %调制=0
Bsim = 1/dt;



y0 = [1.8e18;1e13;1];                                                   %速率方程初值,取初始值 N=Ntr,P=0,此时g=0；相位为1
[t,y] = ode23(@(t,y) rate_equ_nonoise(t,y,I0,Im,f,on,Bsim),tspan, y0);   %不包含噪声 龙格库塔 解微分方程 

I = I0+Im*sin(2*pi*f*t)*on;
%     for j=1:3
%     [t,y] = ode23(@(t,y) rate_equ(t,y,I0,Im,f), tspan, y0);  %计算结果有虚数，换成odde23未能解决
%     y0=[y(length(tspan),1);y(length(tspan),2)];
%     end

% figure(1)
% 
% yyaxis left
% plot(t*1e9,y(:,1))             %ns-N

% % % % % % % % hold on 
% % % % % % % % plot(t*1e9,I*1e3*1e18,'color','red')     %ns-mA*1e18
% % % % % % % % xlim([10 20])
% % % % % % % % hold off
% % % % % % % ylim([0 6])
% % % % % % 
% % % % % % % xlabel('t(ns)','Fontsize',10,'color','black');
% % % % % % % ylabel('N,  Current(\div 10^{18} cm^{-3} mA）','Fontsize',10,'color','black');

% yyaxis right
% plot(t*1e9,y(:,2))
% % ylabel('Photon','Fontsize',10,'color','black');
% 
% title('载流子与光子随时间变化')
% 
% axis square



%% FN谱
% 
% for j=1:1:10        %循环10次取平均
%  L = length(y);
%  d_phi = zeros(L,1);
%  del_v = zeros(L,1);
%  gain = zeros(L,1);
%  for i=1:1:L
%  xi = randn;
%  gain(i) = rate_equ_gain(y(i,1),y(i,2));
%  d_phi(i) = 0.5*alpha*(conf*vg*gain(i)-1/tp)+(1/y(i,2))*sqrt((beta*y(i,2)*B*y(i,1)*y(i,1)*Bsim*0.5))*xi;
%  del_v(i) = d_phi(i)/(2*pi); 
%  end
%  
%  phi = y(:,3);
%  
%  
%  disp(j)
%  
%  start = 3000;          %移除前3000个点（弛豫振荡）
%  N = length(t)- start;  %采样点
%  bin = 0:N; 
%  df=1/(dt);             %采样频率
%  fax_Hz = bin*df/N;     %频率点
%  N_2 = ceil(N/2);

 %=============计算频率噪声===============================================
%  FFT1 = 2.*abs(fft((del_v(start+1:end))))*dt; 
%  FFT1(1)=FFT1(1)/2;
%  FFT1(N_2)=FFT1(N_2)/2;
%  FN = (1/(t2-t(start+1)))*(FFT1.^2);
%  ffp = fax_Hz(1:N_2); 
%  pp = FN(1:N_2);
%  
%  if j==1
%     New = pp;
%     temp = New;
%  else
%     New = (pp + temp)./2;
%     temp = New;
%  end
%  
%  %=============计算相位噪声===============================================
%  
%  FFT2 = 2.*abs(fft((phi(start+1:end))))*dt; %phi ---> del_v
%  FFT2(1)=FFT2(1)/2;
%  FFT2(N_2)=FFT2(N_2)/2;
%  SN = (1/(t2-t(start+1)))*(FFT2.^2);
%  ffs = fax_Hz(1:N_2); 
%  ps = SN(1:N_2);
%  
%  if j==1
%     New1 = ps;
%     temp1 = New1;
%  else
%     New1 = (ps + temp1)./2;
%     temp1 = New1;
%  end
%  
%  
% end
% % 
% % figure(2)
% % plot(t*1e9,del_v);         %频率随时间的变化关系,ns
% % title('频率随时间的变化')
% % xlim([20 30])
% 
% 
% 
% % 
%  figure(3)
%  loglog(ffp, New,'LineWidth',2)
%  title('单边带频率噪声谱')
% ylim([1e4 1e15])
% xlim([1e7 1e10])
% % text(2e7,2e6,'Free Running')
% xlabel('Freqeuncy (Hz)');
% ylabel('Frequency Noise (Hz^2/Hz) ');


%=========================Smooth===========================

% figure(4)
% smoo2 = smooth(ffp, New,0.0001,'loess'); % smoothing function for New
% loglog(ffp, smoo2,'LineWidth',2)
% % title('单边带频率噪声谱')
% ylim([1e5 1e14])
% xlim([2e7 1e10])
% % text(2e7,2e6,'Free Running')
% xlabel('Freqeuncy (Hz)');
% ylabel('Frequency Noise (Hz^2/Hz) ');
% % legend( {'\beta = 10^{-4}'},'FontSize',20)
% hold on
% 
% save Master_Laser_beta_-4
% 
% 
%  
%  figure(5)
%  loglog(ffs, New1 ,'LineWidth',2)
%  title('单边带相位噪声谱')



% %  grid on
% %  set(gcf,'color','w');
% %  xlabel('Freqeuncy (Hz)','Fontsize',40,'color','black');
% %  ylabel('FN (Hz) ','Fontsize',40,'color','black');
% %  set(gca,'Fontsize',40)
% %  Q=get(gca,'ytick');
% %  set(gca,'ytick',Q(2:end))
%%
Nm = y(:,1);
Sm = y(:,2);
Phim = y(:,3);
tm = t;

save Master_nonoise_beta_-6

%% 光谱
mean(Nm)
mean(Sm)

figure 
plot(tm,Sm);

figure
plot(tm,Nm);  

figure
N=length(Sm);               %频域采样点数
Fs=1e12;                        %时域采样频率（Hz），决定了频域的最大值
n=Fs/N;                         %频率域的分辨率 
f0=(-N/2+1:N/2)*n;
Aeff = 2e-4*2e-4;
K = 0.3*vg*Aeff*h*mu;
E = sqrt(Sm).*exp(1i*Phim);
OS = abs( fftshift(fft(ifftshift(sqrt(K)*E))) );
maxa = max(OS);
A = 10*log10(OS./maxa);
plot(f0/1e9,A)
xlim([-30 30])
hold on
ES = abs( fftshift(fft(ifftshift( abs(E).^2 ))) );
maxa = max(ES);
B = 10*log10(ES./maxa);
plot(f0/1e9,B)
