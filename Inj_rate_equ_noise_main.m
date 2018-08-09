clearvars -except Phim Sm tm 



%% 速率方程,光注入+噪声=======对应方程 Inj_rate_equ_noise ======================================
   

%% 激光器参数
q = 1.6e-19;    %C,电子电量
c = 3e10;       %光速度，单位（cm/s）
h = 6.62e-34;  	%普朗克常数
V = 4e-12;      %有源区体积 cm^3
conf=0.032;     %光限制因子
ng=4.2;         %群折射率
vg=c/ng;      	%群速度 cm/s

L = 250*1e-4;   %腔长，cm
width = 2e-4;   %腔宽，cm
d = 80e-8;      %有源区厚度

g0=1800;        %线性增益系数 cm^-1
Ntr=1.8e18;     %透明载流子浓度 cm^-3
Ns=-0.4e18;     %增益参数 cm^-3
eps=1.5e-17;    %增益压缩因子 cm^3

A=0;            %非辐射复合系数 
B=0.8e-10;      %辐射复合系数  cm^3/s
C=3.5e-30;      %俄歇复合系数  cm^6/s
alpha=5;        %线宽增益因子
eta=0.8;        %电流注入效率
beta=0.895e-4;  %自发发射因子

lambda=0.98e-4;	%波长，cm
mu=c/lambda;    %角频率，Hz
tp=2.77e-12;    %光子寿命,s   tau_p 

%% 

t1 = 10e-9;         %Start time  %s微分方程求解时间区间,从td+t1->td+t2
t2 = 50e-9;         %Stop time
dt = 1e-12;         %步长
nsteps = (t2-t1)/(dt) +1;
tspan = linspace(t1,t2,nsteps);    

Im = 2.1e-3;        %调制幅度
f = 2.5e9;            %调制频率，Hz
I0 = 3e-3;          %偏置电流
Bsim = 1/dt;

on = 1;             %调制=0
inj = 1;             %注入锁定=0
%说明 I = I0+Im*sin(2*pi*f*t)*on ;
Kinj = inj*vg/2/L;  % 偏高
scal = 0.001;               %将主激光器的光子数放大，相当于强注入达到稳定的锁定
deta_f = 2e9;         %注入锁定失谐频率设为0

tic
y0 = [1.8e18;1e13;1];                                                                               %速率方程初值,取初始值 N=Ntr,P=0,此时g=0；相位为1
[t,y] = ode23(@(t,y) Inj_rate_equ_noise(t,y,I0,Im,f,on,Bsim,Kinj,Sm,Phim,tm,scal,deta_f),tspan, y0);  	%龙格库塔 解微分方程
toc

% %     for j=1:3
% %     [t,y] = ode23(@(t,y) rate_equ_noise(t,y,I0,Im,f), tspan, y0);  %计算结果有虚数，换成odde23未能解决
% %     y0=[y(length(tspan),1);y(length(tspan),2)];
% %     end

% % figure(1)
% % 
% % yyaxis left
% % plot(t*1e9,y(:,1))             %ns-N
% % % % hold on 
% % % % plot(t*1e9,I*1e3*1e18,'color','red')     %ns-mA*1e18
% % % % xlim([10 20])
% % % % hold off
% % % ylim([0 6])
% % 
% % % xlabel('t(ns)','Fontsize',10,'color','black');
% % % ylabel('N,  Current(\div 10^{18} cm^{-3} mA）','Fontsize',10,'color','black');
% % 
% % yyaxis right
% % plot(t*1e9,y(:,2))
% % % ylabel('Photon','Fontsize',10,'color','black');
% % 
% % title('载流子与光子随时间变化')
% % 
% % axis square



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
%  d_phi(i) = 0.5*alpha*(conf*vg*gain(i)-1/tp) - Kinj*sqrt(scal.*Sm(i)./y(i,2)).*sin(y(i,3)-Phim(i)-2*pi*deta_f*t(i))- 2*pi*deta_f + (1/y(i,2))*sqrt((beta*y(i,2)*B*y(i,1)*y(i,1)*Bsim*0.5))*xi;
%  del_v(i) = d_phi(i)/(2*pi); 
%  end
%  
%  phi = y(:,3);
%  
%  disp(j)
%  
%  start = 3000;          %移除前3000个点（弛豫振荡）
%  N = length(t)- start;  %采样点
%  bin = 0:N; 
%  df=1/(dt);             %采样频率
%  fax_Hz = bin*df/N;     %频率点
%  N_2 = ceil(N/2);
% 
%  %=============计算频率噪声===============================================
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

% % % % figure(2)
% % % % plot(t*1e9,del_v);         %频率随时间的变化关系,ns
% % % % title('频率随时间的变化')
% % % % xlim([20 30])
% % % 
% % % 
% % % % figure(3)
% % % % loglog(ffp, New,'LineWidth',2)
% % % % title('单边带频率噪声谱')




%=========================
% figure(4)
% smoo2 = smooth(ffp, New,0.004,'loess'); % smoothing function for New
% loglog(ffp, smoo2,'LineWidth',2)
% title('单边带频率噪声谱','Fontsize',30)
% set(gca,'FontSize',30);
% ylim([0.2e6 1e12])
% xlim([1e7 1e12])
% % % % text(2e7,2e8,'Free Running')
% xlabel('Freqeuncy (Hz)','Fontsize',30);
% ylabel('Frequency Noise (Hz^2/Hz) ','Fontsize',30);
% hold on
% % % figure(4)
% % % loglog(ffs, New1 ,'LineWidth',2)
% % % title('单边带相位噪声谱')
...

% 
% % %  grid on
% % %  set(gcf,'color','w');
% % %  xlabel('Freqeuncy (Hz)','Fontsize',40,'color','black');
% % %  ylabel('FN (Hz) ','Fontsize',40,'color','black');
% % %  set(gca,'Fontsize',40)
% % %  Q=get(gca,'ytick');
% % %  set(gca,'ytick',Q(2:end))
%%

N_s = y(:,1);
S_s = y(:,2);
Phi_s = y(:,3);
t_s = t;

save sl_mod2.5GHz_inj_detaf2GHz_scal0.001_master_nonoise.mat
%% 光谱
mean(N_s)
mean(S_s)

figure 
plot(t_s,S_s);

figure
plot(t_s,N_s);  

figure
N=length(S_s);               %频域采样点数
Fs=1e12;                        %时域采样频率（Hz），决定了频域的最大值
n=Fs/N;                         %频率域的分辨率 
f0=(-N/2+1:N/2)*n;

Aeff = 2e-4*2e-4;
K = 0.3*vg*Aeff*h*mu;

E = sqrt(S_s).*exp(1i*Phi_s);

OS =abs( fftshift(fft(ifftshift(sqrt(K)*E))) );
maxa = max(OS);
A = 10*log10(OS./maxa);
plot(f0/1e9,A)
xlim([-30 30])
hold on
% 
 ES =abs( fftshift(fft(ifftshift( abs(E).^2 ))) );
 maxa = max(ES);
 B = 10*log10(ES./maxa);
 plot(f0/1e9,B)









