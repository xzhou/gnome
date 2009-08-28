function [] = f(x)
    x(1,1) = 1;
    g(x);
end

function [] = g(x)
    x(2, 2) = 1;
end