function [sx] = copyUpperToLower(x)
%this function copy the upper triangle to lower triangle
[m n] = size(x);
sx = x;
for i = 1:m
    for j = i+1:n
        sx(j,i) = sx(i,j);
    end
end
end