function [tk, err] = OD_LOS_VBR(case_id, if_rand, if_load_LOS)
%% Approach orbit determination using both LOS and range (obtained from the OptNav)
%% Load the Spice ephemeris
cspice_kclear;
load_SPICE_zhou();
%% Parameters setting
load('BennuApproachDataSpice.mat');
if if_load_LOS == 1
    load('Single_OD_Data.mat', 'Euler_real_data');
end
number = 24 * 14;
dt = 3600;
ns = 6;
% Data setting
eleme = zeros(ns, number + 1);
if case_id == 1
    errR = 100;
    errV = 1e-2;
else
    errR = 10;
    errV = 1e-3;
end
if if_rand == 1
    xs0e = xs0 + [randn(3,1) * errR; randn(3,1) * errV];
else
    error00 = [ones(3,1) * errR; ones(3,1) * errV];
    xs0e = xs0 + error00;
end
eleme(:,1) = xs0e;
error0 = [ones(3,1) * errR; ones(3,1) * errV];
P = diag(error0 .* error0);
P0 = P; % the covariance of the estimation at the initial epoch
Q = eye(ns) * 1e-14;
%% Reference orbits
initialEphmerisEpoch = cspice_str2et('2018/9/26 19:15:30');
% [s] % the ephemeris epoch corresponding to initialEpoches(1)
t0 = initialEphmerisEpoch;
tf = number * dt + t0;
%% SPICE propagation
spiceKernelList = {'lsk/naif0011.tls',  'spk/planets/de432s.bsp',  'pck/gm_de431.tpc', ...
    'EM_MoonCenteredRotation.fk',  'EarthCenteredInertial.fk',...
    'SunEarthRotationFrameSCR.fk','SunEarthRotationFrameECR.fk',...
    'MoonCenteredInertial.fk', 'SunCenteredInertial.fk'};
% asteriod dynamic equation in ephemeris model (with SRP)
Aa = 550^2 * pi / 7.3e10 * 1e-6;
fcnSRPa = @(t,X)DynamicEphemerisInertial_SRP(t, X, 'Sun', ...
    {'Sun', 'Mercury', 'Earth', 'Moon', '4'}, spiceKernelList, Aa);
% spacecraft dynamic equation in ephemeris model (with SRP)
As = 8.5 / 2105 * 1e-6;
fcnSRPs = @(t,X)DynamicEphemerisInertial_SRP(t, X, 'Sun', ...
    {'Sun', 'Mercury', 'Earth', 'Moon', '4'}, spiceKernelList, As);
odeOptions = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
%% Reference orbit propagation
tic;
fprintf('================= Asteriod orbit =================\n');
[tk, statea] = ode45(fcnSRPa, [t0: dt: tf], xa0, odeOptions);
time_cost = toc;
fprintf('time cost: %f\n', time_cost);
tic;
fprintf('================= Spacecraft orbit =================\n');
[~, states] = ode45(fcnSRPs, [t0: dt: tf], xs0, odeOptions);
time_cost = toc;
fprintf('time cost: %f\n', time_cost);
obs_data = zeros(3, number + 1); % uesd to save the measurement data
obs_R_data = zeros(3, number + 1);
ratio_data = zeros(1, number + 1);
ratio_fit_data = zeros(1, number + 1);
Aera_data = zeros(1, number + 1);
LOS_data = zeros(3, number + 1);
Euler_data = zeros(2, number + 1);
Center_data = zeros(2, number + 1);
LOS_real_data = zeros(3, number + 1);
if if_load_LOS == 0
    Euler_real_data = zeros(2, number + 1);
end
Euler_err_data = zeros(4, number + 1);
time_cost = zeros(4, number + 1);
cov_data = zeros(ns, number + 1);
Rd_data = zeros(1, number + 1);
x_P_data = struct('xhat', zeros(6,1), 'Phat', zeros(6,6));
for i = 1: number + 1
    x_P_data(i) = struct('xhat', zeros(6,1), 'Phat', zeros(6,6));
end
%% Begin estimation
len = 1;
for i = 1: number + 1
    ti = tk(i) - tk(1);
    fprintf('======== Time: %6f ========\n', ti / 3600 / 24);
    %% Estimated states
    xs = states(i,:);
    xa = statea(i,:);
    xbar = eleme(:,i);
    %% Measurement generation
    image_name = ['Bennu\', num2str(i), '.jpg'];
    Img = imread(image_name);
    [S, ~, ~, center, t1_cost, t2_cost] = CalculateAera_time_cost(Img);
    Aera_data(i) = S;
    Center_data(:,i) = center;
    Euler_rand = 1e-4 * randn(2,1);
    Euler_err_data(:,i) = [center * 13.5e-6; Euler_rand];
    los = (xs(1:3) - xa(1:3)) / norm(xs(1:3) - xa(1:3)); % 标称los
    LOS_data(:,i) = los;
    [az, elev, ~] = cart2sph(los(1), los(2), los(3));
    Euler_data(:,i) = [az; elev];
    if if_load_LOS == 0
        Euler_real_data(:,i) = [az; elev] + center * 13.5e-6 + Euler_rand;
    end
    [los_real_x, los_real_y, los_real_z] = sph2cart(Euler_real_data(1,i), Euler_real_data(2,i), 1);
    LOS_real_data(:,i) = [los_real_x, los_real_y, los_real_z]';
    %% Measurement Update
    tic;
    % Calculate predicted measurement
    if i <= len + 0 % doesn't reverse in the first 10 measurements
        xsei = eleme(:,1);
        xai = statea(1,:)';
        P0i = P0;
        H0 = STM_dis(xsei, xai);
        % Assume that the ephemeris of the asteroid is precise
        [dy, H, R, y, ratio, ratio_fit, Rd] = data_smooth(xs, xa, ...
            xbar, xa, xsei, states(1,:), statea(1,:), Euler_data(:,i), ...
            Euler_real_data(:,i), Aera_data(1), Aera_data(i), ratio_data(:,1: i));
        obs_data(:,i) = y; % save the measurement data
        obs_R_data(:,i) = diag(R);
        ratio_data(:,i) = ratio;
        ratio_fit_data(:,i) = ratio_fit;
        Rd_data(:,i) = Rd;
        H = H(1:2,:);
        dy = dy(1:2,:);
        R = R(1:2,1:2);
    else
        xsei = x_P_data(i-len).xhat;
        P0i = x_P_data(i-len).Phat;
        xai = statea(i-len,:)';
        H0 = STM_dis(xsei, xai);
        d0e = norm(xsei(1:3) - xai(1:3));
        % Assume that the ephemeris of the asteroid is precise
        [dy, H, R, y, ratio, ratio_fit, Rd, de, d0e_] = data_smooth(xs, xa, ...
            xbar, xa, xsei, states(i-len,:), statea(i-len,:), Euler_data(:,i), ...
            Euler_real_data(:,i), Aera_data(i-len), Aera_data(i), ratio_data(:,1: i));
        obs_data(:,i) = y; % save the measurement data
        obs_R_data(:,i) = diag(R);
        ratio_data(:,i) = ratio;
        ratio_fit_data(:,i) = ratio_fit;
        Rd_data(:,i) = Rd;
    end
    %% Calculate the Kalman filter gain
    if length(diag(R)) == 2

    else
        H(end,:) = H(end,:) - H0 * Phi^-1 * de / d0e;
    end
    %% Calculate the estimations
    Pxz = P * H';
    Pzz = H * P * H' + R;
    K_k = Pxz * Pzz^-1;
    xhat = eleme(:,i) + K_k * dy;
    eleme(:,i) = xhat;
    Phat = P - K_k * Pzz * K_k';
    cov_data(:,i) = diag(Phat);
    %% One step prediction
    y0 = [eleme(:,i); reshape(eye(ns), ns^2, 1)];
    [~, yk] = ode45(fcnSRPs, [tk(i), (tk(i) + dt)], y0, odeOptions);
    xbar = yk(end,1: ns)';
    if i <= number
        eleme(:,i + 1) = xbar;
    end
    Phi = reshape(yk(end, ns + 1: end), ns, ns);
    P = Phi * Phat * Phi' + Q;
    x_P_data(i).xhat = xhat;
    x_P_data(i).Phat = Phat;
    t3_cost = toc;
    time_cost(:,i) = [t1_cost; t2_cost; t3_cost; t1_cost + t2_cost + t3_cost];
end
%% Save data
err = eleme(:,1:number+1) - states';
relative_state = states' - statea';
tk = t0: dt: tf;
plot_OD_estimated_error(err, tk, relative_state);