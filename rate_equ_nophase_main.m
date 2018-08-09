clc,clear

%==========���ʷ��̲�����  ���� ��λ  =============================================
%==========��Ӧ���� rate_equ_nophase =======


%% ����

q = 1.6e-19;    %C,���ӵ���
c = 3e10;       %���ٶȣ���λ��cm/s��
h = 6.62e-34; 	%���ʿ˳���

V = 4e-12;      %��Դ����� cm^3
conf = 0.032;     %����������
ng = 4.2;         %Ⱥ������
vg = c/ng;        %Ⱥ�ٶ� cm/s

g0 = 1800;        %��������ϵ��   cm^-1
Ntr = 1.8e18;     %͸��������Ũ�� cm^-3
Ns = -0.4e18;     %�������       cm^-3
eps = 1.5e-17;    %����ѹ������    cm^3

A = 0;            %�Ƿ��临��ϵ�� 
B = 0.8e-10;      %���临��ϵ��  cm^3/s
C = 3.5e-30;      %��Ъ����ϵ��  cm^6/s

alpha = 5;        %�߿���������
eta = 0.8;        %����ע��Ч��
eta_0 = 0.45;     %���ռ�Ч��
beta = 0.895e-4;  %�Է���������

lambda = 0.98e-4;	%������cm
mu = c/lambda;    %��Ƶ�ʣ�Hz
tp = 2.77e-12;    %��������,s
a_i = 5;        %����ģ�cm-1
a_m = 45.6;      %������ģ�cm-1

%% Part1 ֱ��ƫ���£�������,������ ������ı仯,�ɵó� ��ֵ��������ֵ�����ӣ���ֵ����


t1 = 0; 
t2 = 5e-9; 
dt = 5e-14; 
nsteps = (t2-t1)/(dt) +1;
tspan = linspace(t1,t2,nsteps); %span of time values

Im = 1e-3;            %���Ʒ���
f = 2e9;              %����Ƶ�ʣ�Hz
on = 0;               %����=0
I=0:0.1e-3:2e-3;      %����0-20mA
% y0 = [1.8e18;0];
y0 = [10;10];

N_i = zeros(1,length(I));
delta_P_w = zeros(1,length(I));
delta_Pw = zeros(1,length(delta_P_w));
gain_i = zeros(1,length(I));

for i=1:1:length(I)
    I0 = I(i);
    for j=1:3
    [t,y] = ode23(@(t,y) rate_equ_nophase(t,y,I0,Im,f,on), tspan, y0);  %������������������odde23δ�ܽ��
    y0=[y(length(tspan),1);y(length(tspan),2)];
    end
    L=length(y); 
    N_i(i) = y(L,1);          %��ͬ��������̬������
    delta_P_w(i) = y(L,2);    %��ͬ��������̬������
end

delta_Pw = eta_0*h*mu.*delta_P_w*(V/conf)*vg*a_m;         	%����
gain_i = g0./(1+eps.*delta_P_w).*log((N_i+Ns)/(Ntr+Ns));      %����

figure                        %������,������ ������ı仯

yyaxis left
plot(1000*I,N_i/(1e18),'-o')
text(0.6,1.5,'g')
hold on 
plot(1000*I,1000*delta_Pw*10,'-d','color','red')
text(0.5,3.2,'N')
xlabel('Current(mA)','Fontsize',10,'color','black');
ylabel('N(\times 10^{18} cm^{-3}),  P(\times 0.1mW��','Fontsize',10,'color','black');

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
% plot(N_i,gain_i,'-o'); 	%�����������ӵı仯
% xlabel('Carrier density (cm^{-3})','Fontsize',10,'color','black');
% ylabel('gain(cm^{-1})','Fontsize',10,'color','black');
% grid on
% axis square


%% Part2 ֱ��ƫ�û��ߵ���״̬�� �������������ʱ��仯��ϵ

t1 = 10e-9;     	%Start time  %s΢�ַ������ʱ������,��td+t1->td+t2
t2 = 100e-9;        %Stop time
dt = 5e-12;         %����
nsteps = (t2-t1)/(dt) +1;
tspan = linspace(t1,t2,nsteps);  


f = 2e9;            %����Ƶ�ʣ�Hz
I0 = 3e-3;          %ƫ�õ���
Im = 2e-3;          %���Ʒ���
on = 0             %����=0



y0 = [1.8e18;1e13];                                                 
[t,y] = ode23(@(t,y) rate_equ(t,y,I0,Im,f,on),tspan, y0);   
I = I0+Im*sin(2*pi*f*t)*on;
    for j=1:3
    [t,y] = ode23(@(t,y) rate_equ_nophase(t,y,I0,Im,f), tspan, y0);  %������������������odde23δ�ܽ��
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

%% Part3 ��������С�źŵ��Ʋ���

%ֱ��ƫ��,�ض������»�ȡ��̬������,С�źŲ���================================


% f = 2e9;            %����Ƶ�ʣ�Hz
% I0 = 3e-3;          %ƫ�õ���
% Im = I0-1e-3;       %���Ʒ���
% 
% 
% on = 0;                       %����=0
% 
% t1 = 1e-9; % Start time       %΢�ַ������ʱ������,��td+t1->td+t2
% t2 = 5e-9; % Stop time
% dt = 5e-14;                   %����
% nsteps = (t2-t1)/(dt) +1;
% tspan = linspace(t1,t2,nsteps); %span of time values
% y0 = [1e18;0];
% 
% % for j=1:3
% % [t,y] = ode45(@(t,y) rate_equ(t,y,I0,Im,f,on), tspan, y0);  %������������������odde23δ�ܽ��
% % y0=[y(length(tspan),1);y(length(tspan),2)];
% % end
% [t,y] = ode45(@(t,y) rate_equ(t,y,I0,Im,f,on), tspan, y0);
% I = I0+Im*sin(2*pi*f*t)*on;
% 
% L=length(y);
% N = y(L,1)                    %��̬������
% P = y(L,2)                    %��̬������
% figure
% plot(t,y(:,2))                %�������ȶ�����
% 
% %����N,Pֵ����С�źŲ���===================================================
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
% f = 0:0.01e9:10e9;            %0-10GHz������Ӧ
% H = zeros(1,length(f));
% 
% H = wr^2./(wr^2-(2*pi*f).^2+1i*2*pi*f*r);
% Hw = 10*log10(abs(H).^2);
% 
% figure
% plot(f/1e9,Hw)                %С�ź���Ӧ���� ������GHz
% xlabel('Frequency(GHz)','Fontsize',10,'color','black');
% ylabel('Respones','Fontsize',10,'color','black');
% grid on
% axis square


%% Part4 ���ʷ��̽�С�ź���Ӧ

% % % I0 = 3e-3;%ƫ�õ���
% % % Im = 0.5e-3;%���Ʒ���
% I = I0+Im*sin(2*pi*f*t)*on;

% % % on = 1;                       %����=1
% % % 
% % % 
% % % mod_f=0.2e9:0.1e9:10e9;      %����Ƶ��
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
%     f_max = find(localmaxmin(ph,'max')); %����ֵ��λ��
%     jmax = ph(f_max);
%     f_min = find(localmaxmin(ph,'min'));
%     jmin = ph(f_min);
% 
%     delta_P_w(i) = 0.5*(  mean(jmax) - mean(jmin) );       %��ͬ����Ƶ���£��������仯�ķ���
% end


% Hw = delta_P_w/Im/(eta*conf*tp/q/V);                    %��ͬ����Ƶ���£�������Ӧ,���й�ʽ��5.56��
% R = 20*log10(abs(Hw));                                  %��ͬ����Ƶ���£�������Ӧ
% 
% figure(1)               %���Ƶ������ߣ���Ӧ�Ĺ������仯���ߣ���ʱ����Ƶ��Ϊ max(mode_f)

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


% figure(2)               %������Ӧ����
% plot(mod_f/1e9,R,'-o');     
% xlabel('Frequency(GHz)','Fontsize',10,'color','black');
% ylabel('Response','Fontsize',10,'color','black');
% grid on
% axis square
% plot(mod_f/1e9,delta_P_w,'-o'); 

%%


















