%% 主程序rate_equ_nophase_main 计算特定电流下的稳态值后，计算小信号响应

N %稳态载流子数
P %稳态光子数
gain = rate_equ_gain(N,P);
a0 = 5.3386e-16;  %cm^2
ap = eps*gain/(1+eps*P);
a = a0/(1+eps*P);

rnn = A+2*B*N+3*C*N^2+P*vg*a;
rnp = 1/conf/tp-vg*ap*P;
rpn = conf*vg*a*P;
rpp = conf*vg*ap*P;

wr = (rnp*rpn+rnn*rpp)^0.5;
fr = wr/(2*pi)  %Hz
r = rnn+rpp     %Hz

f = 0:0.01e9:10e9;
H = zeros(1,length(w));

H = wr^2./(wr^2-(2*pi*f).^2+1i*2*pi*f*r);
Hw = log10(abs(H).^2);

plot(f/1e9,Hw)  %% 横坐标GHz

axis square