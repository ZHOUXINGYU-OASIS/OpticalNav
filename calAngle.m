function theta = calAngle(r1, r2)
%% �����������н�
cos_theta = r1' * r2 / norm(r1) / norm(r2);
theta = acos(cos_theta) * 180 / pi;