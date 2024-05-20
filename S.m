function [outputArg1] = S(inputArg1)
%INTERPOLATION3 此处显示有关此函数的摘要
%   此处显示详细说明
w = abs(inputArg1);
if w<1
    S = 1-2*w^2+w^3;
elseif (1<=w) && (w<=2)
    S = -8*w+5*w^2+w^3;
else
    S = 0;
end
outputArg1 = S;