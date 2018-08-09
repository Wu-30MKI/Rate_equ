function gain=rate_equ_gain(N,P)
% g0=1284;          %��������ϵ��cm^-1
% Ntr=1.2e18;       %͸��������Ũ��cm^-3
% Ns=0.92*Ntr;      %�������cm^-3
% eps=18/Ntr;       %����ѹ������cm^3

g0=1800;            %��������ϵ��cm^-1
Ntr=1.8e18;         %͸��������Ũ��cm^-3
Ns=-0.4e18;         %�������cm^-3
eps=1.5e-17;        %����ѹ������cm^3
gain=g0/(1+eps*P)*log((N+Ns)/(Ntr+Ns));
