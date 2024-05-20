function [S, Igray, Imax, center, t1_cost, t2_cost] = CalculateAera_time_cost(I)
%% Calculate the area also record the time cost
[M00, M11R, M11I, M20] = Zernike77();
if length(size(I))==3
    I = rgb2gray(I); 
end
Imax = max(max(I));
%% Calculate the center
I = fig01(I);
[x, y] = find(I);
x_max = max(x) + 5;
x_min = min(x) - 5;
y_max = max(y) + 5;
y_min = min(y) - 5;
I = I(x_min: x_max, y_min: y_max);
tic;
[X, Y] = center_Moment_My(I);
center = [X + x_min - 540.5; Y + y_min - 540.5]; % 1080P
t1_cost = toc;
tic;
Igray = graythresh(I);
I = im2bw(I, 0.1);
K = double(I);
% Convolution operation
A11I = conv2(M11I, K);
A11R = conv2(M11R, K);
A20 = conv2(M20, K);
A00 = conv2(M00, K);
% Cut off the excess
A11I = A11I(4:end-3,4:end-3);
A11R = A11R(4:end-3,4:end-3);
A20 = A20(4:end-3,4:end-3);
A00 = A00(4:end-3,4:end-3);
J = zeros(size(K));
theta = atan2(A11I, A11R); % calculate theta
% Calculate the three parameters of the edge
A11C = A11R .* cos(theta) + A11I .* sin(theta);
l = A20 ./ A11C;
k = 1.5 * A11C ./ ((1 - l.^2) .^ 1.5);
e = abs(l) > 1 / 3.5;
k(e) = 0;
g = graythresh(k);
% Edge judgment condition
a = abs(l) < 1 / sqrt(2) * 2 / 7;
b = abs(k) > g;
% a and b are the judging results of distance and edge strength respectively
J(a & b) = 1;
I41 = imfill(J,'holes');
I4 = bwperim(I41);
J = bwmorph(I4,'thin',Inf);
[x, y] = find(J == 1); % Pixel level coordinates of the edge
Z = [x + l(find(J == 1)) .* cos(theta(find(J == 1))), ...
     y + l(find(J == 1)) .* sin(theta(find(J == 1)))];
% Subpixel coordinate
x = Z(:,1);
y = Z(:,2);
[x1, y1] = Fun_ReSort(x, y);
S = polyarea(x1,y1);
t2_cost = toc;