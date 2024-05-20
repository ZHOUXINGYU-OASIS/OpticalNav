function [pred, delta] = ratio_polynomial_smooth(x, y, n)
%% smooth the ratio data using the convex optimization
[nx, mx] = size(x);
l = nx;
if nx < mx
    x = x';
    l = mx;
end
[ny, my] = size(y);
if ny < my
    y = y';
end
A = zeros(l, n + 1);
for i = 1: n + 1
    A(:,i) = x.^(i - 1);
end
Q = [A,  -eye(l, l);
     -A, -eye(l, l);
     A, zeros(l, l)];
k_vector = [y; -y; ones(l,1)];
c_vector = [zeros(n + 1, 1); ones(l, 1)];
param.MSK_IPAR_LOG = 0;
res = msklpopt(c_vector, Q, -ones(length(k_vector), 1) * inf, k_vector, ...
    -ones(length(c_vector), 1) * inf, [inf; 0; ones(length(c_vector) - 2, 1) * inf], param, 'minimize');
THETA = res.sol.itr.xx;
pred = A(end,:) * THETA(1: (n+1));
pred_all = A(:,:) * THETA(1: (n+1));
normr = norm(pred_all - y);
df = max(0, length(y) - (n+1));
%% predict the std
R = triu(A);
V = A(end,:);
E = V / R;
e = sqrt(1 + sum(E.*E,2));
delta = std(pred_all - y);