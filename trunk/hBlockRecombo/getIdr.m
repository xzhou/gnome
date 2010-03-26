function [idr] = getIdr(caseSeq, testSeq, fdr, caseR, refR)
%% get the identification rate
TrCase = getTrM(caseSeq, caseR, refR);
TrTest = getTrM(testSeq, caseR, refR);
T = getThreshold(TrTest, fdr);
idr = sum(TrCase > T);
end