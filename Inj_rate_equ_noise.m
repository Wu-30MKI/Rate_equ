
%% ���ʷ��̣���ע������+����================================================


%% 

function dy = Inj_rate_equ_noise(t,y,I0,Im,f,on,Bsim,Kinj,Sm,Phim,tm,scal,deta_f)

q = 1.6e-19;    %C,���ӵ���
c = 3e10;       %���ٶȣ���λ��cm/s��
h = 6.62e-34; 	%���ʿ˳���
V = 4e-12;      %��Դ����� cm^3
conf=0.032;     %����������
ng=4.2;         %Ⱥ������
vg=c/ng;        %Ⱥ�ٶ� cm/s

g0=1800;        %��������ϵ��   cm^-1
Ntr=1.8e18;     %͸��������Ũ�� cm^-3
Ns=-0.4e18;     %�������       cm^-3
eps=1.5e-17;    %����ѹ������    cm^3

A=0;            %�Ƿ��临��ϵ�� 
B=0.8e-10;     	%���临��ϵ��  cm^3/s
C=3.5e-30;      %��Ъ����ϵ��  cm^6/s
alpha=5;        %�߿���������
eta=0.8;        %����ע��Ч��
beta=0.895e-4;      %�Է���������

lambda=0.98e-4;	%������cm
mu=c/lambda;     %��Ƶ�ʣ�Hz
tp=2.77e-12;     %��������,s   tau_p

%% �������

xis=randn;
xin=randn;
xi=randn;

Dnn = sqrt((2*(A*y(1)+B*y(1)*y(1)+C*y(1)*y(1)*y(1))*Bsim))*(xin);
Dss = sqrt((beta*2*y(2)*B*y(1)*y(1)*Bsim))*(xis);
Dpp = sqrt((beta*y(2)*B*y(1)*y(1)*Bsim*0.5))*(xi);

Sml = interp1(tm, Sm, t);       % ��ֵ
Phiml = interp1(tm, Phim, t);   % ��ֵ

I=I0+Im*sin(2*pi*f*t)*on;
%% ���ʷ��� 

g=rate_equ_gain(y(1),y(2));
dy=zeros(3,1);
dy(1) = eta*I/q/V-(A*y(1)+B*(y(1))^2+C*(y(1))^3)-vg*g*y(2)+Dnn-Dss;
dy(2) = (conf*vg*g-1/tp)*y(2)+conf*beta*B*(y(1))^2+2*Kinj*sqrt(scal*Sml.*y(2)).*cos(y(3)-Phiml-2*pi*deta_f*t)+Dss;
dy(3) = alpha/2*(conf*vg*g-1/tp)-Kinj*sqrt(scal*Sml./y(2)).*sin(y(3)-Phiml-2*pi*deta_f*t)  +(1/y(2))*Dpp;
              % ������֣���� dy(3) ��ʾ����slave��λ������Ӧ��ȥ�� - 2*pi*deta_f
              % ���շ��̵��㷨������Master field �ı�ʾΪ Em = E*exp��1i*2*pi*delta_f*t��
%%









