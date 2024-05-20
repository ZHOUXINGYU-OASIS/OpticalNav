function H = cal_H_angle_dis_ratio(xs, xa)
%% 计算测角+测速度距离比值的H矩阵
rsx = xs(1);
rsy = xs(2);
rsz = xs(3);
vsx = xs(4);
vsy = xs(5);
vsz = xs(6);
rax = xa(1);
ray = xa(2);
raz = xa(3);
vax = xa(4);
vay = xa(5);
vaz = xa(6);
H = zeros(4, 6);
% lx
dlx_drsx = 1/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2) - (abs(rax - rsx)*sign(rax - rsx)*(rax - rsx))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
dlx_drsy = -(abs(ray - rsy)*sign(ray - rsy)*(rax - rsx))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
dlx_drsz = -(abs(raz - rsz)*sign(raz - rsz)*(rax - rsx))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
% ly
dly_drsx = -(abs(rax - rsx)*sign(rax - rsx)*(ray - rsy))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
dly_drsy = 1/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2) - (abs(ray - rsy)*sign(ray - rsy)*(ray - rsy))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
dly_drsz = -(abs(raz - rsz)*sign(raz - rsz)*(ray - rsy))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
% lz
dlz_drsx = -(abs(rax - rsx)*sign(rax - rsx)*(raz - rsz))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
dlz_drsy = -(abs(ray - rsy)*sign(ray - rsy)*(raz - rsz))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
dlz_drsz = 1/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2) - (abs(raz - rsz)*sign(raz - rsz)*(raz - rsz))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);

dr_drsx = (2*abs(rax - rsx)*sign(rax - rsx)*((conj(vax) - conj(vsx))*(rax - rsx) + (conj(vay) - conj(vsy))*(ray - rsy) + (conj(vaz) - conj(vsz))*(raz - rsz)))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^2 - (conj(vax) - conj(vsx))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2);
dr_drsy = (2*abs(ray - rsy)*sign(ray - rsy)*((conj(vax) - conj(vsx))*(rax - rsx) + (conj(vay) - conj(vsy))*(ray - rsy) + (conj(vaz) - conj(vsz))*(raz - rsz)))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^2 - (conj(vay) - conj(vsy))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2);
dr_drsz = (2*abs(raz - rsz)*sign(raz - rsz)*((conj(vax) - conj(vsx))*(rax - rsx) + (conj(vay) - conj(vsy))*(ray - rsy) + (conj(vaz) - conj(vsz))*(raz - rsz)))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^2 - (conj(vaz) - conj(vsz))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2);
dr_dvsx = -(rax - rsx)/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2);
dr_dvsy = -(ray - rsy)/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2);
dr_dvsz = -(raz - rsz)/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2);
H(1:3, 1:3) = [dlx_drsx, dlx_drsy, dlx_drsz;
             dly_drsx, dly_drsy, dly_drsz;
             dlz_drsx, dlz_drsy, dlz_drsz];
H(4, 1:6) = [dr_drsx, dr_drsy, dr_drsz, dr_dvsx, dr_dvsy, dr_dvsz];