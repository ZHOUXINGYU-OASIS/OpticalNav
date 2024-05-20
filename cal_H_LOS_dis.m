function H = cal_H_LOS_dis(xs, xa)
%% Calculate the Jacobi matrix of the angle and distance measurements
rsx = xs(1);
rsy = xs(2);
rsz = xs(3);
rax = xa(1);
ray = xa(2);
raz = xa(3);
H = zeros(4, 6);
%% Angle measurement
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
%% Distance measurement
% dis
dr_drsx = -(abs(rax - rsx)*sign(rax - rsx))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2);
dr_drsy = -(abs(ray - rsy)*sign(ray - rsy))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2); 
dr_drsz = -(abs(raz - rsz)*sign(raz - rsz))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2);
%% H matrix
H(1:3, 1:3) = [dlx_drsx, dlx_drsy, dlx_drsz;
               dly_drsx, dly_drsy, dly_drsz;
               dlz_drsx, dlz_drsy, dlz_drsz];
H(4, 1:3) = [dr_drsx, dr_drsy, dr_drsz];