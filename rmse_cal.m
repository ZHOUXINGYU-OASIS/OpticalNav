function r = rmse_cal(x)
r = sqrt(sum(x.^2) / length(x));