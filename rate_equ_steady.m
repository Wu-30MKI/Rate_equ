clc,clear

%% 用解析公式算激光器参数 Ith gth,Nth,wr,r

%% 参数

q = 1.6e-19;        %C,电子电量
c = 3e10;           %光速度，单位（cm/s）
h = 6.62e-34;       %普朗克常数
V = 4e-12;          %有源区体积 cm^3
conf=0.032;         %光限制因子
ng=4.2;             %群折射率
vg=c/ng;            %群速度 cm/s

g0=1800;            %线性增益系数   cm^-1
Ntr=1.8e18;         %透明载流子浓度 cm^-3
Ns=-0.4e18;         %增益参数       cm^-3
eps=1.5e-17;        %增益压缩因子    cm^3

A=0;                %非辐射复合系数 
B=0.8e-10;          %辐射复合系数  cm^3/s
C=3.5e-30;          %俄歇复合系数  cm^6/s
alpha=5;            %线宽增益因子
eta=0.8;            %电流注入效率
beta=0.869e-4;      %自发发射因子

lambda=0.98e-4;     %波长，cm
mu=c/lambda;        %角频率，Hz
tp=2.77e-12;        %光子寿命,s

a_i=5;              %内损耗 cm-1
a_m=1/vg/tp-a_i;    %cm-1
%% L-I 

Nth=Ntr*exp(1/(conf*vg*tp*g0))
Ith=q*V*(A*Nth+B*Nth^2+C*Nth^3)/eta
gth=1/conf/vg/tp
I=0:0.01e-3:2e-3;
P=zeros(1,length(I));
for i=1:1:length(I)
    if I(i)<=Ith
            P(i)=B*Nth^2/(A*Nth+B*Nth^2+C*Nth^3)*eta*a_m/(a_m+a_i)*h*mu/q*beta*I(i);
    else
            P(i)=eta*a_m/(a_m+a_i)*h*mu/q*(I(i)-Ith);
   
    end
end
plot(1000*I,1000*P)  %mW-mA 
%% 

