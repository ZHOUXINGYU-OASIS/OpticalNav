function d_SPR_dr = diff_SPR_dx(r, rs)
% 所有输入都是坐标转换后的状态
x = r(1);
y = r(2);
z = r(3);
xs = rs(1);
ys = rs(2);
zs = rs(3);
SPR_x_rx = 1/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (3*abs(x - xs)*sign(x - xs)*(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2);
SPR_x_ry = -(3*abs(y - ys)*sign(y - ys)*(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2);
SPR_x_rz = -(3*abs(z - zs)*sign(z - zs)*(x - xs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2);
SPR_y_rx = -(3*abs(x - xs)*sign(x - xs)*(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2);
SPR_y_ry = 1/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (3*abs(y - ys)*sign(y - ys)*(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2);
SPR_y_rz = -(3*abs(z - zs)*sign(z - zs)*(y - ys))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2);
SPR_z_rx = -(3*abs(x - xs)*sign(x - xs)*(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2);
SPR_z_ry = -(3*abs(y - ys)*sign(y - ys)*(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2);
SPR_z_rz = 1/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(3/2) - (3*abs(z - zs)*sign(z - zs)*(z - zs))/(abs(x - xs)^2 + abs(y - ys)^2 + abs(z - zs)^2)^(5/2);
d_SPR_dr = [SPR_x_rx, SPR_x_ry, SPR_x_rz;
            SPR_y_rx, SPR_y_ry, SPR_y_rz;
            SPR_z_rx, SPR_z_ry, SPR_z_rz];