function gain=rate_equ_gain(N,P)
% g0=1284;          %线性增益系数cm^-1
% Ntr=1.2e18;       %透明载流子浓度cm^-3
% Ns=0.92*Ntr;      %增益参数cm^-3
% eps=18/Ntr;       %增益压缩因子cm^3

g0=1800;            %线性增益系数cm^-1
Ntr=1.8e18;         %透明载流子浓度cm^-3
Ns=-0.4e18;         %增益参数cm^-3
eps=1.5e-17;        %增益压缩因子cm^3
gain=g0/(1+eps*P)*log((N+Ns)/(Ntr+Ns));
