function [ lc00 ] = matrix2c00( M )
%MATRIX2C00 Summary of this function goes here
%   Detailed explanation goes here
[m n] = size(M);
c00 = ones(n, n);
lc00 = [];
for i = 1:n-1
    for j = (i+1):n
        snp1 = M(:,i);
        snp2 = M(:,j);
        ac00 = sum((snp1 + snp2) == 2);
        c00(i,j) = ac00;
        %linearIndex = calcPairIndex(i,j,n);
        lc00 = [lc00; ac00];
    end
end
end
