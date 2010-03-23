function [estR] = hapSeqEstR(hapSeq01)
%calculate estimate R from haplotype sequence, 1 = major, 0 = minor
genoSeq = combineHapSeq(hapSeq01);
estR = estimateR(genoSeq);
end
