function [maskMatrix] = getMaskP(n, p)
%create a n*n mask with 1-p percent elements equals to -1
maskMatrix = ones(n);
for k = 1:n-1
    x = (rand(n-k, 1)<p)*2 -1;
    maskMatrix(k,k+1:n) = x;
end
maskMatrix = copyUpperToLower(maskMatrix);
end