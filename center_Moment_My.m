function [X, Y] = center_Moment_My(g)
%% 灰度中心的质心提取
% Independent Optical Navigation Processing for the Osiris-Rex 
% Mission Using the Goddard Image Analysis and Navigation Tool
g = double(g);
[n, m] = size(g);
xc = [0; 0];
for i = 1: n
    for j = 1: m
        xc = xc + [i - 1; j - 1] * g(i,j);
%         xc = xc + [j; i] * g(i,j);
    end
end
xc = xc / sum(sum(g));
X = xc(1);
Y = xc(2);