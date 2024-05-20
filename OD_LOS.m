function [tk, err] = OD_LOS(case_id, if_rand, if_load_LOS)
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
%% Data setting
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
Q = eye(ns) * 1e-14;
covData = zeros(ns, number + 1);
covData(:,1) = diag(P);
%% Initial state
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
tic;
fprintf('================= Asteriod orbit =================\n');
[~, state_asteriod] = ode45(fcnSRPa, [t0: dt: tf], xa0, odeOptions);
time_cost = toc;
fprintf('time cost: %f\n', time_cost);
tic;
fprintf('================= Spacecraft orbit =================\n');
[~, elem] = ode45(fcnSRPs, [t0: dt: tf], xs0, odeOptions);
time_cost = toc;
fprintf('time cost: %f\n', time_cost);
%% Begin estimation
i = 1;
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
for t = t0: dt: tf
    fprintf('time:%f\n', (t - t0) / 3600 / 24);
    %% Measurement generation
    xs = elem(i,:)';
    xa = state_asteriod(i,:)';
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
    tic;
    %% Observation
    [dy, H, R, y, h] = EKF_angle_true(elem(i,:)', state_asteriod(i,:)', ...
        eleme(:,i), state_asteriod(i,:)', 1e-4, Euler_real_data(:,i));
    %% 计算矩阵增益
    K_k = P * H' * (H * P * H' + R)^-1;
    %% 计算估计值
    eleme(:,i) = eleme(:,i) + K_k * dy;
    P = (eye(ns) - K_k * H) * P;
    covData(:,i) = diag(P);
    y0 = [eleme(:,i); reshape(eye(ns), ns^2, 1)];
    [~, yk] = ode45(fcnSRPs, [t, (t + dt)], y0, odeOptions);
    %% One step prediction
    if i <= number
        eleme(:,i + 1) = yk(end,1: ns)';
        Phi = reshape(yk(end, ns + 1: end), ns, ns);
        P = Phi * P * Phi' + Q;
    end
    t3_cost = toc;
    time_cost(:,i) = [t1_cost; t2_cost; t3_cost; t1_cost + t2_cost + t3_cost];
    i = i + 1;
end
err = eleme - elem';
relative_state = elem' - state_asteriod';
tk = t0: dt: tf;
%% Plot figure
plot_OD_estimated_error(err, tk - t0, relative_state);