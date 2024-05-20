function fig = fig01_k(fig, k)
[n, m] = size(fig);
for i = 1: n
    for j = 1: m
        if fig(i,j) <= k
            fig(i,j) = 0;
        end
    end
end