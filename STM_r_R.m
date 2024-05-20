function [H, HR] = STM_r_R(xs, xa, R)
%% 对于小行星图像尺寸角的测量STM矩阵
rsx = xs(1);
rsy = xs(2);
rsz = xs(3);
rax = xa(1);
ray = xa(2);
raz = xa(3);
H = zeros(1, 6);
dalpha_drsx = (abs(rax - rsx)*sign(rax - rsx)*real(R))/((real(R)^2 + (imag(R) - (abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2))^2)*(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2));
dalpha_drsy = (abs(ray - rsy)*sign(ray - rsy)*real(R))/((real(R)^2 + (imag(R) - (abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2))^2)*(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2));
dalpha_drsz = (abs(raz - rsz)*sign(raz - rsz)*real(R))/((real(R)^2 + (imag(R) - (abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2))^2)*(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2));
H(1:3) = [dalpha_drsx, dalpha_drsy, dalpha_drsz];
HR = -(imag(R) - (abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2))/(real(R)^2 + (imag(R) - (abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2))^2);