function [ErrData, TimeData, stdData, rmseData, time_cost, MDData] = MC_OD_LOS_VBR5(case_id)
%% Load the Spice ephemeris
cspice_kclear;
load_SPICE_zhou();
%% Parameters setting
Us = 1.32712440e+11;
load('BennuApproachDataSpice.mat');
number = 24 * 14;
dt = 3600;
ns = 6;
%% Begin MC
Num = 100;
ErrData = zeros(Num, ns * 2, number + 1);
TimeData = zeros(Num, 4, number + 1);
MDData = zeros(Num, number + 1);
dn = 5;
for i_MC = 1: Num
    fprintf('MC:%d\n', i_MC);
    %$ Data setting
    eleme = zeros(ns, number + 1);
    if case_id == 2
        errR = 10;
        errV = 1e-3;
    end
    if case_id == 1
        errR = 100;
        errV = 1e-2;
    end
    xs0e = xs0 + [randn(3,1) * errR; randn(3,1) * errV];
    eleme(:,1) = xs0e;
    error0 = [ones(3,1) * errR; ones(3,1) * errV];
    P = diag(error0 .* error0);
    P0 = P; % the covariance of the estimation at the initial epoch
    Q = eye(ns) * 1e-14;
    MDdata = zeros(1, number + 1);
    MDdata(1, 1) = cal_MD(xs0, xs0e, P);
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
    Euler_real_data = zeros(2, number + 1);
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
        %% Estimated states
        xs = states(i,:);
        xa = statea(i,:);
        %% Measurement generation
        image_name = ['Bennu\', num2str(i), '.jpg'];
        Img = imread(image_name);
        [S, Igray, Imax, center, t1_cost, t2_cost] = CalculateAera_time_cost(Img);
        Aera_data(i) = S;
        Center_data(:,i) = center;
        Euler_rand = 1e-4 * randn(2,1);
        Euler_err_data(:,i) = [center * 13.5e-6; Euler_rand];
        los = (xs(1:3) - xa(1:3)) / norm(xs(1:3) - xa(1:3)); % 标称los
        LOS_data(:,i) = los;
        [az, elev, ~] = cart2sph(los(1), los(2), los(3));
        Euler_data(:,i) = [az; elev];
        Euler_real_data(:,i) = [az; elev] + center * 13.5e-6 + Euler_rand;
        [los_real_x, los_real_y, los_real_z] = sph2cart(Euler_real_data(1,i), Euler_real_data(2,i), 1);
        LOS_real_data(:,i) = [los_real_x, los_real_y, los_real_z]';
        %% Measurement Update
        tic;
        % Calculate predicted measurement
        if i <= len + 0 % doesn't reverse in the first 10 measurements
            % Assume that the ephemeris of the asteroid is precise
            [dy, H, R, y, ratioData] ...
                = data_smooth_multi(xs, xa, ...
                eleme(:,1:i), statea(1:i,:)', x_P_data(1:i), Euler_data(:,i), ...
                Euler_real_data(:,i), Aera_data(1:i), 1, []);
        else
            % Assume that the ephemeris of the asteroid is precise
            y0_reverse = [eleme(:,i); reshape(eye(ns), ns^2, 1)];
            if i <= dn
                l_s = 1:i;
                tk_reverse = tk(l_s);
                tk_reverse = tk_reverse(end:-1:1);
                [~, yk_reverse] = ode45(fcnSRPs, tk_reverse, y0_reverse, odeOptions);
                if length(tk_reverse) == 2
                    yk_reverse = yk_reverse([1,end],:);
                end
                yk_reverse = yk_reverse(end:-1:1,:);
                [dy, H, R, y, ratioData] ...
                    = data_smooth_multi(xs, xa, ...
                    eleme(:,l_s), statea(l_s,:)', x_P_data(l_s), Euler_data(:,i), ...
                    Euler_real_data(:,i), Aera_data(1:i), length(l_s), yk_reverse(:,:));
            else
                l_s = (i - dn):i;
                tk_reverse = tk(l_s);
                tk_reverse = tk_reverse(end:-1:1);
                [~, yk_reverse] = ode45(fcnSRPs, tk_reverse, y0_reverse, odeOptions);
                if length(tk_reverse) == 2
                    yk_reverse = yk_reverse([1,end],:);
                end
                yk_reverse = yk_reverse(end:-1:1,:);
                [dy, H, R, y, ratioData] ...
                    = SREIF_angle_multi_dis_measurement_prediction_xs0_smooth(xs, xa, ...
                    eleme(:,l_s), statea(l_s,:)', x_P_data(l_s), Euler_data(:,i), ...
                    Euler_real_data(:,i), Aera_data(1:i), length(l_s), yk_reverse(:,:));
            end
        end
        %% Calculate the Kalman filter gain
        K_k = P * H' * (H * P * H' + R)^-1;
        %% Calculate the estimations
        xhat = eleme(:,i) + K_k * dy;
        eleme(:,i) = xhat;
        Phat = (eye(ns) - K_k * H) * P;
        cov_data(:,i) = diag(Phat);
        %% One step prediction
        y0 = [eleme(:,i); reshape(eye(ns), ns^2, 1)];
        [~, yk] = ode45(fcnSRPs, [tk(i), (tk(i) + dt)], y0, odeOptions);
        xbar = yk(end,1: ns)';
        if i <= number
            MDdata(1, i + 1) = cal_MD(eleme(:,i), states(i, :)', P);
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
    len = length(tk);
    err_RF = err;
    for i = 1: len
        a = eye(3);
        a(:, 1) = relative_state(1:3, i);
        b = Schmidt_orthogonalization(a);
        err_RF(1:3, i) = b' * err(1:3, i);
        err_RF(4:6, i) = b' * err(4:6, i);
    end
    ErrData(i_MC,:,:) = [err; err_RF];
    TimeData(i_MC,:,:) = time_cost;
    MDData(i_MC, :) = MDdata;
end
stdData = zeros(number + 1, ns * 2);
for i = 1: number + 1
    for j = 1: ns * 2
        stdData(i,j) = std(ErrData(:,j,i));
    end
end
rmseData = zeros(number + 1, ns);
for i = 1: number + 1
    for j = 1: ns * 2
        rmseData(i,j) = rmse_cal(ErrData(:,j,i));
    end
end
time_cost = mean(mean(TimeData, 3), 1);