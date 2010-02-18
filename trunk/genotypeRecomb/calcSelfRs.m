function [rs, r] = calcSelfRs(doubleHapSeq, alleleMapping)
%the input data is two haplotyoe sequence, we calculate the "rsquare"
[r c00 c01 c10 c11] = calcPairwiseFreq(doubleHapSeq, alleleMapping);
[m n] = size(doubleHapSeq);
r = zeros(n, n);
rs = r;

for i = 1:n
    for j = i:n
        r(i,j) = (c00(i,j) + c11(i,j) - c01(i,j) - c10(i,j)+2)/4;
        r(j,i) = r(i,j);
    end
end
rs = r.*r;
end

