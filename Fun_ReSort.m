function [xnew,ynew] = Fun_ReSort(x,y)
%% [xnew,ynew] = Fun_ReSort(x,y) ������ɢ������x-y�����վ������ԭ������Ϊ��������xnew-ynew
%   ��ʼ��Ϊx��y�ĵ�1��x(1)-y(1)

% ��ʼ��
xnew = zeros(length(x), 1);
ynew = zeros(length(y), 1);

p = [x(1), y(1)];       % ��ǰ��
Ind = 2:length(x);     % ʣ����±� 
ii = 1;
xnew(ii) = p(1);
ynew(ii) = p(2);

while length(Ind) > 1
    [~,ind] = min((p(1) - x(Ind)).^2 + (p(2) - y(Ind)).^2);
    p = [x(Ind(ind)) y(Ind(ind))];    % ��ǰ��
    Ind = setdiff(Ind,Ind(ind));      % ʣ����±� 
    ii = ii+1;
    xnew(ii) = p(1);
    ynew(ii) = p(2);
end

xnew(end) = x(Ind);  % ���һ��
ynew(end) = y(Ind);