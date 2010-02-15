function [newMatrix] = setDiag(x, val)
newMatrix = x;
[m n] = size(x);
if m ~= n
    newMatrix = x;
end
for i = 1:m
    newMatrix(i,i) = val;
end
end