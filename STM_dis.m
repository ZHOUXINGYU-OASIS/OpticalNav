function H = STM_dis(xs, xa)
%% 计算仅测距的H矩阵
rsx = xs(1);
rsy = xs(2);
rsz = xs(3);
rax = xa(1);
ray = xa(2);
raz = xa(3);
H = zeros(1, 6);
dd_drsx = -(abs(rax - rsx)*sign(rax - rsx))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2);
dd_drsy = -(abs(ray - rsy)*sign(ray - rsy))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2);
dd_drsz = -(abs(raz - rsz)*sign(raz - rsz))/(abs(rax - rsx)^2 + abs(ray - rsy)^2 + abs(raz - rsz)^2)^(1/2);
H(:, 1:3) = [dd_drsx, dd_drsy, dd_drsz];