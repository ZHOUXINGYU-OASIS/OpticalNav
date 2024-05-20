function [dy, H, R, y, h] = EKF_angle(xs, xa, xse, xae, obs_std)
%% calculate the angle measurement in EKF
[l, ~] = measurement_angle_dis(xs, xa);
y = l + randn(3,1) * obs_std;
[le, ~] = measurement_angle_dis(xse, xae);
h = le;
dy = y - h;
%% Jacobi matrix
H = cal_H_angle_dis(xse, xae);
H = H(1:3,:);
%% Measurement covariance
R = diag(ones(1,3) * obs_std^2);