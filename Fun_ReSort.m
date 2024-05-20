function [xnew,ynew] = Fun_ReSort(x,y)
%% [xnew,ynew] = Fun_ReSort(x,y) 将无序散点序列x-y，按照距离最近原则重排为有序序列xnew-ynew
%   起始点为x和y的第1点x(1)-y(1)

% 初始化
xnew = zeros(length(x), 1);
ynew = zeros(length(y), 1);

p = [x(1), y(1)];       % 当前点
Ind = 2:length(x);     % 剩余点下标 
ii = 1;
xnew(ii) = p(1);
ynew(ii) = p(2);

while length(Ind) > 1
    [~,ind] = min((p(1) - x(Ind)).^2 + (p(2) - y(Ind)).^2);
    p = [x(Ind(ind)) y(Ind(ind))];    % 当前点
    Ind = setdiff(Ind,Ind(ind));      % 剩余点下标 
    ii = ii+1;
    xnew(ii) = p(1);
    ynew(ii) = p(2);
end

xnew(end) = x(Ind);  % 最后一点
ynew(end) = y(Ind);