function [l, d] = measurement_angle_dis(xs, xa)
%% calculate the measurement vector
rsx = xs(1);
rsy = xs(2);
rsz = xs(3);
rax = xa(1);
ray = xa(2);
raz = xa(3);
rs = [rsx; rsy; rsz];
ra = [rax; ray; raz];
d = norm(rs - ra);
lx = (rsx - rax) / d;
ly = (rsy - ray) / d;
lz = (rsz - raz) / d;
l = [lx; ly; lz];
[az, elev, ~] = cart2sph(l(1), l(2), l(3));
l = [az; elev];