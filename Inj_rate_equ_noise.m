
%% 速率方程，光注入锁定+噪声================================================


%% 

function dy = Inj_rate_equ_noise(t,y,I0,Im,f,on,Bsim,Kinj,Sm,Phim,tm,scal,deta_f)

q = 1.6e-19;    %C,电子电量
c = 3e10;       %光速度，单位（cm/s）
h = 6.62e-34; 	%普朗克常数
V = 4e-12;      %有源区体积 cm^3
conf=0.032;     %光限制因子
ng=4.2;         %群折射率
vg=c/ng;        %群速度 cm/s

g0=1800;        %线性增益系数   cm^-1
Ntr=1.8e18;     %透明载流子浓度 cm^-3
Ns=-0.4e18;     %增益参数       cm^-3
eps=1.5e-17;    %增益压缩因子    cm^3

A=0;            %非辐射复合系数 
B=0.8e-10;     	%辐射复合系数  cm^3/s
C=3.5e-30;      %俄歇复合系数  cm^6/s
alpha=5;        %线宽增益因子
eta=0.8;        %电流注入效率
beta=0.895e-4;      %自发发射因子

lambda=0.98e-4;	%波长，cm
mu=c/lambda;     %角频率，Hz
tp=2.77e-12;     %光子寿命,s   tau_p

%% 随机噪声

xis=randn;
xin=randn;
xi=randn;

Dnn = sqrt((2*(A*y(1)+B*y(1)*y(1)+C*y(1)*y(1)*y(1))*Bsim))*(xin);
Dss = sqrt((beta*2*y(2)*B*y(1)*y(1)*Bsim))*(xis);
Dpp = sqrt((beta*y(2)*B*y(1)*y(1)*Bsim*0.5))*(xi);

Sml = interp1(tm, Sm, t);       % 插值
Phiml = interp1(tm, Phim, t);   % 插值

I=I0+Im*sin(2*pi*f*t)*on;
%% 速率方程 

g=rate_equ_gain(y(1),y(2));
dy=zeros(3,1);
dy(1) = eta*I/q/V-(A*y(1)+B*(y(1))^2+C*(y(1))^3)-vg*g*y(2)+Dnn-Dss;
dy(2) = (conf*vg*g-1/tp)*y(2)+conf*beta*B*(y(1))^2+2*Kinj*sqrt(scal*Sml.*y(2)).*cos(y(3)-Phiml-2*pi*deta_f*t)+Dss;
dy(3) = alpha/2*(conf*vg*g-1/tp)-Kinj*sqrt(scal*Sml./y(2)).*sin(y(3)-Phiml-2*pi*deta_f*t)  +(1/y(2))*Dpp;
              % 推算后发现，如果 dy(3) 表示的是slave相位，所以应该去掉 - 2*pi*deta_f
              % 按照方程的算法，这里Master field 的表示为 Em = E*exp（1i*2*pi*delta_f*t）
%%









