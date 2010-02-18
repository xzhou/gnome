function [f] = singleAlleleFreqOfHapSeq(hapSeq4Encode, majorAllele)
% get the single allele frequency of each snps
hapSeq01 = zeroOneEncodeHapSeq(hapSeq4Encode, majorAllele);
[m n] = size(hapSeq01);
f = (m - sum(hapSeq01))/m;
f = f';
end

