function [dy, H, R, y, ratioData] = data_smooth_multi(xs, xa, xse_data, xae_data, P_data, los, los_real, S_data, len, yk_reverse)
%% The function to compute the measurement angle and distance
%% Using the true LOS and range measurements
%% Using ratio data smoothing
%% Using multi-distance data
% dy = zeros(1 + len, 1);
h = zeros(1 + len, 1); % predicted measurements
y = zeros(1 + len, 1); % income measurements (with noise)
H = zeros(1 + len, 6); % Jacobi matrix
R = zeros(1 + len, 1 + len); % covariance matrix of measurement noises
len_k = 10;
l_S = length(S_data);
xse = xse_data(:,end);
xae = xae_data(:,end);
ratioData = zeros(len - 1, 1);
for k = 1: len - 1
    xse0 = xse_data(:,end-k);
    xse0_reverse = yk_reverse(end-k,1:6)';
    Phi_reverse = reshape(yk_reverse(end-k,7:42), 6, 6);
    xae0 = xae_data(:,end-k);
    P0 = P_data(end-k).Phat;
    H0 = STM_dis(xse0_reverse, xae0);
    d0e = norm(xse0_reverse(1:3)' - xae0(1:3)'); % the estimated intial distance
    if l_S >= k + len_k
        ratio_data = zeros(len_k,1);
        for i = 1: len_k
            ratio_data(i) = sqrt(S_data(end + 1 - i - k) / S_data(end + 1 - i));
        end
        ratio_data = ratio_data(end:-1:1);
        id = 1:length(ratio_data);
        [ratio_fit, err] = ratio_polynomial_smooth(id, ratio_data(id), 1);
    else
        ratio_data = zeros(l_S - k,1);
        for i = 1: (l_S - k)
            ratio_data(i) = sqrt(S_data(end + 1 - i - k) / S_data(end + 1 - i));
        end
        ratio_data = ratio_data(end:-1:1);
        if length(ratio_data) <= 5
            ratio_fit = ratio_data(end);
            err = 0.008;
%             err = 0.05;
        else
            id = 1:length(ratio_data);
            [ratio_fit, err] = ratio_polynomial_smooth(id, ratio_data(id), 1);
        end
    end
    ratioData(k) = ratio_fit;
    de = d0e * ratio_fit;
    y(k) = de;
    h(k) = norm(xse(1:3)' - xae(1:3)'); % the estimated intial distance
    R(k,k) = (d0e * err)^2;
    Hk = STM_dis(xse, xae) - H0 * Phi_reverse * h(k) / d0e;
    H(k,:) = Hk;
end
y(end-1:end) = los_real;
[le, ~] = measurement_angle_dis(xse, xae);
h(end-1:end) = le;
Hk = STM_angle(xse, xae);
H(end-1:end,:) = Hk(1:2,:);
R(end-1:end,end-1:end) = diag([1.2e-4; 1.03e-4].^2);
dy = y - h;
if dy(end-1) > pi
    dy(end-1) = dy(end-1) - 2 * pi;
end
if dy(end-1) < -pi
    dy(end-1) = dy(end-1) + 2 * pi;
end