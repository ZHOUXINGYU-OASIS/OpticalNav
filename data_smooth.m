function [dy, H, R, y, ratio_e, ratio_fit, Rd, de, d0e] = data_smooth(xs, xa, xse, xae, xse0, xs0, xa0, los, los_real, S0, S, ratio_data)
%% The function to compute the measurement angle and distance
%% Using the true LOS and range measurements
%% Using ratio data smoothing
ratio_e = sqrt(S0 / S);
ratio_fit = ratio_e;
ratio_data(end) = ratio_e;
err = ratio_e * 0.005;
% err = ratio_e * 0.005;
if length(ratio_data) > 10
    %% smoothing
    id = (length(ratio_data) - 10) : length(ratio_data);
    [ratio_fit, err] = ratio_polynomial_smooth(id, ratio_data(id), 1);
elseif length(ratio_data) <= 10 && length(ratio_data) >= 5
    %% smoothing
    id = 1 : length(ratio_data);
    [ratio_fit, err] = ratio_polynomial_smooth(id, ratio_data(id), 1);
end
%% True measurements
d0e = norm(xse0(1:3)' - xa0(1:3)); % the estimated intial distance
d = norm(xs(1:3) - xa(1:3)); % the true distance at the current epoch
[~, d_real] = measurement_angle_dis_d0(xs, xa, d0e, ratio_fit);
l = los_real;
y = [l; d_real];
y_ref = [los; d]; % reference measurements
%% Predicted measurements
[le, de] = measurement_angle_dis(xse, xae);
h = [le; de];
dy = y - h;
if dy(1) > pi
    dy(1) = dy(1) - 2 * pi;
end
if dy(1) < -pi
    dy(1) = dy(1) + 2 * pi;
end
%% Jacobi matrix
H = [STM_angle(xse, xae); STM_dis(xse, xae)];
%% Measurement covariance
% R = diag([y_ref(1:2) - y(1:2); d0e * err].^2);
% R = diag([y_ref(1:3) - y(1:3)].^2);
R = diag([1.2e-4; 1.03e-4; d0e * err * 1].^2);
% R = diag([1.03e-4; 1.2e-4; d0e * 0.0085].^2);
Rd = (y_ref(3) - y(3))^2;