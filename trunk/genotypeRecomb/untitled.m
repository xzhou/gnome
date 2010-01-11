x = zeros(10,10);
k = 0;

for i = 1:9
    for j = i+1:10
        x(i,j) = k;
        k = k + 1;
    end
end

for i = 1:9
    x(i+1:10, i) = x(i, i+1:10);
end

x

