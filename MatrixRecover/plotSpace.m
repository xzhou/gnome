function plotSpace()
    close all;
    for L = 8
        figure;
        ratio = [];
        N = 2:50;
        for i = 1:size(N,2)
            ratio = [ratio, spaceRatioExact(N(i), L)];
        end
        line(N, ratio);
    end
end