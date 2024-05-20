function d_SPR_dr = diff_SPR_dr(r, rm, rs)
% 所有输入都是坐标转换后的状态
rx = r(1);
ry = r(2);
rz = r(3);
rmx = rm(1);
rmy = rm(2);
rmz = rm(3);
rsx = rs(1);
rsy = rs(2);
rsz = rs(3);
SPR_x_rx = 709261484669871/(104857600*(abs(rsx - rmx + rx)^2 + abs(rsy - rmy + ry)^2 + abs(rsz - rmz + rz)^2)^(3/2)) - (3*abs(rsx - rmx + rx)*sign(rsx - rmx + rx)*((3546307423349355*rsx)/524288 - (3546307423349355*rmx)/524288 + (3546307423349355*rx)/524288))/(1000*(abs(rsx - rmx + rx)^2 + abs(rsy - rmy + ry)^2 + abs(rsz - rmz + rz)^2)^(5/2));
SPR_x_ry = -(3*abs(rsy - rmy + ry)*sign(rsy - rmy + ry)*((3546307423349355*rsx)/524288 - (3546307423349355*rmx)/524288 + (3546307423349355*rx)/524288))/(1000*(abs(rsx - rmx + rx)^2 + abs(rsy - rmy + ry)^2 + abs(rsz - rmz + rz)^2)^(5/2));
SPR_x_rz = -(3*abs(rsz - rmz + rz)*sign(rsz - rmz + rz)*((3546307423349355*rsx)/524288 - (3546307423349355*rmx)/524288 + (3546307423349355*rx)/524288))/(1000*(abs(rsx - rmx + rx)^2 + abs(rsy - rmy + ry)^2 + abs(rsz - rmz + rz)^2)^(5/2));
SPR_y_rx = -(3*abs(rsx - rmx + rx)*sign(rsx - rmx + rx)*((3546307423349355*rsy)/524288 - (3546307423349355*rmy)/524288 + (3546307423349355*ry)/524288))/(1000*(abs(rsx - rmx + rx)^2 + abs(rsy - rmy + ry)^2 + abs(rsz - rmz + rz)^2)^(5/2));
SPR_y_ry = 709261484669871/(104857600*(abs(rsx - rmx + rx)^2 + abs(rsy - rmy + ry)^2 + abs(rsz - rmz + rz)^2)^(3/2)) - (3*abs(rsy - rmy + ry)*sign(rsy - rmy + ry)*((3546307423349355*rsy)/524288 - (3546307423349355*rmy)/524288 + (3546307423349355*ry)/524288))/(1000*(abs(rsx - rmx + rx)^2 + abs(rsy - rmy + ry)^2 + abs(rsz - rmz + rz)^2)^(5/2));
SPR_y_rz = -(3*abs(rsz - rmz + rz)*sign(rsz - rmz + rz)*((3546307423349355*rsy)/524288 - (3546307423349355*rmy)/524288 + (3546307423349355*ry)/524288))/(1000*(abs(rsx - rmx + rx)^2 + abs(rsy - rmy + ry)^2 + abs(rsz - rmz + rz)^2)^(5/2));
SPR_z_rx = -(3*abs(rsx - rmx + rx)*sign(rsx - rmx + rx)*((3546307423349355*rsz)/524288 - (3546307423349355*rmz)/524288 + (3546307423349355*rz)/524288))/(1000*(abs(rsx - rmx + rx)^2 + abs(rsy - rmy + ry)^2 + abs(rsz - rmz + rz)^2)^(5/2));
SPR_z_ry = -(3*abs(rsy - rmy + ry)*sign(rsy - rmy + ry)*((3546307423349355*rsz)/524288 - (3546307423349355*rmz)/524288 + (3546307423349355*rz)/524288))/(1000*(abs(rsx - rmx + rx)^2 + abs(rsy - rmy + ry)^2 + abs(rsz - rmz + rz)^2)^(5/2));
SPR_z_rz = 709261484669871/(104857600*(abs(rsx - rmx + rx)^2 + abs(rsy - rmy + ry)^2 + abs(rsz - rmz + rz)^2)^(3/2)) - (3*abs(rsz - rmz + rz)*sign(rsz - rmz + rz)*((3546307423349355*rsz)/524288 - (3546307423349355*rmz)/524288 + (3546307423349355*rz)/524288))/(1000*(abs(rsx - rmx + rx)^2 + abs(rsy - rmy + ry)^2 + abs(rsz - rmz + rz)^2)^(5/2));
d_SPR_dr = [SPR_x_rx, SPR_x_ry, SPR_x_rz;
            SPR_y_rx, SPR_y_ry, SPR_y_rz;
            SPR_z_rx, SPR_z_ry, SPR_z_rz];