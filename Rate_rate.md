速率方程函数：
​	

**rate_equ_nophase :** 	

​			最简单的速率方程，光子密度与载流子密度的关系，没有相位项

​			用来仿真静态和动态特性

​			对应主程序 rate_equ_nophase_main

**rate_equ_nonoise :**	

​			包括相位项的速率方程，不包括噪声项

​			用来产生的光注入时的 主激光器参数 Sm Nm Phim (没有考虑噪声，所以主激光器更加接近理想，还没有算)

​			对应主程序 rate_equ_nonoise_main

**rate_equ_noise   :**	

​			包括噪声的速率方程，当然包含相位项

​			用来产生的光注入时的 主激光器参数 Sm Nm Phim（自发发射因子小对应噪声低）

​			用来计算没有光注入时的 【 光谱，频谱，相噪，频率噪声，RIN 】

​			对应主程序 rate_equ_noise_main

**Inj_rate_equ_noise :**	

​			光注入速率方程，包括噪声	

​			用来计算光注入后的 相噪，频率噪声，RIN

​			对应主程序 Inj_rate_equ_noise_main

