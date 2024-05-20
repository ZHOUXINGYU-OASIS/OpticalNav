function dy = dynamics2BP(t, y, mu)
%% Two-body dynamics
r = y(1:3);
v = y(4:6);
R = norm(r);
dy(1) = v(1);
dy(2) = v(2);
dy(3) = v(3);
dy(4) = -mu / R^3 * r(1);
dy(5) = -mu / R^3 * r(2);
dy(6) = -mu / R^3 * r(3);
dy = dy';