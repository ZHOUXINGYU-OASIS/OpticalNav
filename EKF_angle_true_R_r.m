function [dy, H, R, y, h] = EKF_angle_true_R_r(xse, xae, y, Re)
%% calculate the angle measurement in EKF
[le, de] = measurement_angle_dis(xse, xae);
R_r_e = atan2(Re, de);
h = [le; R_r_e];
dy = y - h;
%% Jacobi matrix
[H1, H2] = STM_r_R(xse, xae, Re);
H = [STM_angle(xse, xae); H1];
H = [H, [0; 0; H2]];
%% Measurement covariance
R = diag(([1.2, 1.03, 13.5e-2] * 1e-4) .^ 2);