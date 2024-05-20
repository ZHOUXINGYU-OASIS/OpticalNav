function H = STM_angle(xs, xa)
%% 计算仅测角的H矩阵
rsx = xs(1);
rsy = xs(2);
rsz = xs(3);
rax = xa(1);
ray = xa(2);
raz = xa(3);
H = zeros(2, 6);
daz_drsx = (imag(rax) - imag(rsx) + real(ray) - real(rsy))/((imag(rax) - imag(rsx) + real(ray) - real(rsy))^2 + (imag(ray) - imag(rsy) - real(rax) + real(rsx))^2);
daz_drsy = (imag(ray) - imag(rsy) - real(rax) + real(rsx))/((imag(rax) - imag(rsx) + real(ray) - real(rsy))^2 + (imag(ray) - imag(rsy) - real(rax) + real(rsx))^2);
daz_drsz = 0;
delev_drsx = -((imag((2*rax - 2*rsx)/((rax - rsx)^2 + (ray - rsy)^2)^(1/2))/(2*(real(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) + imag(raz) - imag(rsz))) - (real((2*rax - 2*rsx)/((rax - rsx)^2 + (ray - rsy)^2)^(1/2))*(imag(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) - real(raz) + real(rsz)))/(2*(real(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) + imag(raz) - imag(rsz))^2))*(real(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) + imag(raz) - imag(rsz))^2)/((real(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) + imag(raz) - imag(rsz))^2 + (imag(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) - real(raz) + real(rsz))^2);
delev_drsy = -((imag((2*ray - 2*rsy)/((rax - rsx)^2 + (ray - rsy)^2)^(1/2))/(2*(real(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) + imag(raz) - imag(rsz))) - (real((2*ray - 2*rsy)/((rax - rsx)^2 + (ray - rsy)^2)^(1/2))*(imag(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) - real(raz) + real(rsz)))/(2*(real(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) + imag(raz) - imag(rsz))^2))*(real(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) + imag(raz) - imag(rsz))^2)/((real(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) + imag(raz) - imag(rsz))^2 + (imag(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) - real(raz) + real(rsz))^2);
delev_drsz = (real(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) + imag(raz) - imag(rsz))/((real(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) + imag(raz) - imag(rsz))^2 + (imag(((rax - rsx)^2 + (ray - rsy)^2)^(1/2)) - real(raz) + real(rsz))^2);
H(:, 1:3) = [daz_drsx, daz_drsy, daz_drsz;
             delev_drsx, delev_drsy, delev_drsz];