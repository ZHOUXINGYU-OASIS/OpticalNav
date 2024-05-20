function dfx = DFx(r, mu)
%% 对天体引力求偏导数
%       x' = f(x)
%       dfx = df / dx
x = r(1);
y = r(2);
z = r(3);
normr = norm(r);
fxx = mu * (3 * x^2 - normr^2) / normr^5;
fyy = mu * (3 * y^2 - normr^2) / normr^5;
fzz = mu * (3 * z^2 - normr^2) / normr^5;
fxy = 3 * mu * x * y / normr^5;
fyx = 3 * mu * y * x / normr^5;
fxz = 3 * mu * x * z / normr^5;
fzx = 3 * mu * z * x / normr^5;
fyz = 3 * mu * y * z / normr^5;
fzy = 3 * mu * z * y / normr^5;
dfx = [fxx fxy fxz;
       fyx fyy fyz;
       fzx fzy fzz];