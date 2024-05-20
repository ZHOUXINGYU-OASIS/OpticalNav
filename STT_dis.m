function H = STT_dis(xs, xa)
%% 计算仅测距的H矩阵
rsx = xs(1);
rsy = xs(2);
rsz = xs(3);
rax = xa(1);
ray = xa(2);
raz = xa(3);
H = zeros(1, 6, 6);
dd_drsx_drsx = sign(rax - rsx)^2/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2) - (abs(rax - rsx)^2*sign(rax - rsx)^2)/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2) + (2*abs(rax - rsx)*dirac(rax - rsx))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2);
dd_drsx_drsy = -(abs(rax - rsx)*abs(ray - rsy)*sign(rax - rsx)*sign(ray - rsy))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
dd_drsx_drsz = -(abs(rax - rsx)*abs(raz - rsz)*sign(rax - rsx)*sign(raz - rsz))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
dd_drsy_drsx = -(abs(rax - rsx)*abs(ray - rsy)*sign(rax - rsx)*sign(ray - rsy))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
dd_drsy_drsy = sign(ray - rsy)^2/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2) - (abs(ray - rsy)^2*sign(ray - rsy)^2)/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2) + (2*abs(ray - rsy)*dirac(ray - rsy))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2);
dd_drsy_drsz = -(abs(ray - rsy)*abs(raz - rsz)*sign(ray - rsy)*sign(raz - rsz))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
dd_drsz_drsx = -(abs(rax - rsx)*abs(raz - rsz)*sign(rax - rsx)*sign(raz - rsz))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
dd_drsz_drsy = -(abs(ray - rsy)*abs(raz - rsz)*sign(ray - rsy)*sign(raz - rsz))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2);
dd_drsz_drsz = sign(raz - rsz)^2/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2) - (abs(raz - rsz)^2*sign(raz - rsz)^2)/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(3/2) + (2*abs(raz - rsz)*dirac(raz - rsz))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2);
H(1, 1:3, 1:3) = [dd_drsx_drsx, dd_drsx_drsy, dd_drsx_drsz;
                  dd_drsy_drsx, dd_drsy_drsy, dd_drsy_drsz;
                  dd_drsz_drsx, dd_drsz_drsy, dd_drsz_drsz];