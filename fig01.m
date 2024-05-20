function fig = fig01(fig)
[n, m] = size(fig);
for i = 1: n
    for j = 1: m
        if fig(i,j) <= 5
            fig(i,j) = 0;
        end
    end
end