function plotSpace()
    N = 10;
    ratio = [];
    L = 2:20;
    for i = 1:size(L,2)
        ratio = [ratio, spaceRatioExact(N, L(i))];
    end
    plot(L, ratio);
    
end