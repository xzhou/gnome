function [seq2int] = zeroOneEncodeHapSeq(seq4int, majorAllele)
%this function will convert the 0~3 encoding of haplotype sequence to 01
%encoding after knowing major allele mapping
[m n] = size(seq4int);
seq2int = (seq4int == repmat(majorAllele, m, 1)) + 0;
end