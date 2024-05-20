function [dy, H, R, y, h] = EKF_angle_true(xs, xa, xse, xae, obs_std, y)
%% calculate the angle measurement in EKF
% [l, ~] = measurement_angle_dis(xs, xa);
% y = l + randn(3,1) * obs_std;
[le, ~] = measurement_angle_dis(xse, xae);
h = le;
dy = y - h;
%% Jacobi matrix
H = STM_angle(xse, xae);
% H = cal_H_angle_dis(xse, xae);
H = H(1:2,:);
%% Measurement covariance
R = diag(([1.2, 1.03] * 1e-4) .^ 2);
% R = diag(([1.03, 1.2] * 1e-4) .^ 2);