function dy = stateTrans(t, y, mu)
%% STM
n = 6; % r1 v1 also need to be propagated
X = y(1: n);
Phi = reshape(y(n + 1: end), n, n);
dX = dynamics2BP(t, X, mu);
%% STM'
dfX = DFx(X(1:3), mu); % r2 v2
dfX = [zeros(3, 3), eye(3, 3);
       dfX,         zeros(3, 3)];
dPhi = dfX * Phi;
dy = [dX; reshape(dPhi, n * n, 1)];