clc,clear

%% �ý�����ʽ�㼤�������� Ith gth,Nth,wr,r

%% ����

q = 1.6e-19;        %C,���ӵ���
c = 3e10;           %���ٶȣ���λ��cm/s��
h = 6.62e-34;       %���ʿ˳���
V = 4e-12;          %��Դ����� cm^3
conf=0.032;         %����������
ng=4.2;             %Ⱥ������
vg=c/ng;            %Ⱥ�ٶ� cm/s

g0=1800;            %��������ϵ��   cm^-1
Ntr=1.8e18;         %͸��������Ũ�� cm^-3
Ns=-0.4e18;         %�������       cm^-3
eps=1.5e-17;        %����ѹ������    cm^3

A=0;                %�Ƿ��临��ϵ�� 
B=0.8e-10;          %���临��ϵ��  cm^3/s
C=3.5e-30;          %��Ъ����ϵ��  cm^6/s
alpha=5;            %�߿���������
eta=0.8;            %����ע��Ч��
beta=0.869e-4;      %�Է���������

lambda=0.98e-4;     %������cm
mu=c/lambda;        %��Ƶ�ʣ�Hz
tp=2.77e-12;        %��������,s

a_i=5;              %����� cm-1
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

