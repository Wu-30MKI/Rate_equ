%% ���ʷ���
%====== ���������������������ע��,������λ===��Ӧ���� rate_equ_nonoise==============
%====== ���� �����ע���µ� ��̬���� �� ��̬���� ==========================
%====== ͨ���ı� �Է���������beta���õ����������Ĳ��� Sm Nm Phim ======================
% clc 
% clear
% clf
% close all
%% ����������
q = 1.6e-19;    %C,���ӵ���
c = 3e10;       %���ٶȣ���λ��cm/s��
h = 6.62e-34;  	%���ʿ˳���
V = 4e-12;      %��Դ����� cm^3
conf=0.032;     %����������
ng=4.2;         %Ⱥ������
vg=c/ng;      	%Ⱥ�ٶ� cm/s

g0=1800;        %��������ϵ�� cm^-1
Ntr=1.8e18;     %͸��������Ũ�� cm^-3
Ns=-0.4e18;     %������� cm^-3
eps=1.5e-17;    %����ѹ������ cm^3

A=0;            %�Ƿ��临��ϵ�� 
B=0.8e-10;      %���临��ϵ��  cm^3/s
C=3.5e-30;      %��Ъ����ϵ��  cm^6/s
alpha=5;        %�߿���������
eta=0.8;        %����ע��Ч��
beta=0.895e-6;      %�Է���������

lambda=0.98e-4;	%������cm
mu=c/lambda;    %��Ƶ�ʣ�Hz
tp=2.77e-12;    %��������,s   tau_p

%% ֱ��ƫ���£�������,������ ������ı仯,�ɵó� ��ֵ��������ֵ�����ӣ���ֵ����

t1 = 0; 
t2 = 5e-9; 
dt = 5e-14;
nsteps = (t2-t1)/(dt) +1;
tspan = linspace(t1,t2,nsteps); 
Bsim = 1/dt;

Im = 1e-3;            %���Ʒ���
f = 2e9;              %����Ƶ�ʣ�Hz
on = 0;               %����=0
I=0:0.1e-3:1e-3;      %����0-10mA
% y0 = [1.8e18;1e13;1]; 
y0 = [Ntr;1e15;1];

N_i = zeros(1,length(I));           %����������
P_i = zeros(1,length(I));           %��������
Power = zeros(1,length(P_i));       %���ʣ�W
gain_i = zeros(1,length(I));

for i=1:1:length(I)
    I0 = I(i);
    [t,y] = ode23(@(t,y) rate_equ_nonoise(t,y,I0,Im,f,on,Bsim),tspan, y0);
%     for j=1:3                                            
%     [t,y] = ode23(@(t,y) rate_equ_switch(t,y,I0,Im,f,on,Bsim),tspan, y0);
%     y0=[y(length(tspan),1);y(length(tspan),2);y(length(tspan))];
%     end
    L=length(y); 
    N_i(i) = y(L,1);          %��ͬ��������̬������
    P_i(i) = y(L,2);    %��ͬ��������̬������
end

Power = eta_0*h*mu.*P_i*(V/conf)*vg*a_m;               %����
% % gain_i = g0./(1+eps.*delta_P_w).*log((N_i+Ns)/(Ntr+Ns));        %����


figure            %ֱ��ƫ���£�������,������������ı仯

yyaxis left
plot(1000*I,N_i/(1e18),'-o')
text(0.6,1.5,'g')
hold on 
plot(1000*I,1000*Power*10,'-d','color','red')
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



% figure(2)                 %�����������ӵı仯
% plot(N_i,gain_i,'-o'); 	
% xlabel('Carrier density (m^(-3))','Fontsize',10,'color','black');
% ylabel('gain(cm^{-1})','Fontsize',10,'color','black');
% grid on
% axis square



%% ֱ��ƫ�û��ߵ���״̬�� �������������ʱ��仯��ϵ

t1 = 10e-9;     	%Start time  %s΢�ַ������ʱ������,��td+t1->td+t2
t2 = 50e-9;        %Stop time
dt = 1e-12;         %����
nsteps = (t2-t1)/(dt) +1;
tspan = linspace(t1,t2,nsteps);    

Im = 1e-3;       %���Ʒ���
f = 5e9;            %����Ƶ�ʣ�Hz
I0 = 3e-3;          %ƫ�õ���
on = 0;             %����=0
Bsim = 1/dt;



y0 = [1.8e18;1e13;1];                                                   %���ʷ��̳�ֵ,ȡ��ʼֵ N=Ntr,P=0,��ʱg=0����λΪ1
[t,y] = ode23(@(t,y) rate_equ_nonoise(t,y,I0,Im,f,on,Bsim),tspan, y0);   %���������� ������� ��΢�ַ��� 

I = I0+Im*sin(2*pi*f*t)*on;
%     for j=1:3
%     [t,y] = ode23(@(t,y) rate_equ(t,y,I0,Im,f), tspan, y0);  %������������������odde23δ�ܽ��
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
% % % % % % % ylabel('N,  Current(\div 10^{18} cm^{-3} mA��','Fontsize',10,'color','black');

% yyaxis right
% plot(t*1e9,y(:,2))
% % ylabel('Photon','Fontsize',10,'color','black');
% 
% title('�������������ʱ��仯')
% 
% axis square



%% FN��
% 
% for j=1:1:10        %ѭ��10��ȡƽ��
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
%  start = 3000;          %�Ƴ�ǰ3000���㣨��ԥ�񵴣�
%  N = length(t)- start;  %������
%  bin = 0:N; 
%  df=1/(dt);             %����Ƶ��
%  fax_Hz = bin*df/N;     %Ƶ�ʵ�
%  N_2 = ceil(N/2);

 %=============����Ƶ������===============================================
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
%  %=============������λ����===============================================
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
% % plot(t*1e9,del_v);         %Ƶ����ʱ��ı仯��ϵ,ns
% % title('Ƶ����ʱ��ı仯')
% % xlim([20 30])
% 
% 
% 
% % 
%  figure(3)
%  loglog(ffp, New,'LineWidth',2)
%  title('���ߴ�Ƶ��������')
% ylim([1e4 1e15])
% xlim([1e7 1e10])
% % text(2e7,2e6,'Free Running')
% xlabel('Freqeuncy (Hz)');
% ylabel('Frequency Noise (Hz^2/Hz) ');


%=========================Smooth===========================

% figure(4)
% smoo2 = smooth(ffp, New,0.0001,'loess'); % smoothing function for New
% loglog(ffp, smoo2,'LineWidth',2)
% % title('���ߴ�Ƶ��������')
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
%  title('���ߴ���λ������')



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

%% ����
mean(Nm)
mean(Sm)

figure 
plot(tm,Sm);

figure
plot(tm,Nm);  

figure
N=length(Sm);               %Ƶ���������
Fs=1e12;                        %ʱ�����Ƶ�ʣ�Hz����������Ƶ������ֵ
n=Fs/N;                         %Ƶ����ķֱ��� 
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
