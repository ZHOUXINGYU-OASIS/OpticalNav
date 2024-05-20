function MD = cal_MD(xr, xe, P)
error = xr - xe;
MD = sqrt(error' * P^(-1) * error);