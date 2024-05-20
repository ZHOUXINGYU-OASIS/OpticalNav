function dy = solar_dynamics(t,y)
%% 航天器在太阳系中飞行的动力学

Us = 1.32712440e+11;

r = y(1:3);
v = y(4:6);

R = norm(r);

dy(1) = v(1);
dy(2) = v(2);
dy(3) = v(3);
dy(4) = -Us/R^3*r(1);
dy(5) = -Us/R^3*r(2);
dy(6) = -Us/R^3*r(3);

dy = dy';