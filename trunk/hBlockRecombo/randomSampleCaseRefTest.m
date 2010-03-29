function [refSeq, caseSeq, testSeq] = randomSampleCaseRefTest(hapSeq, nCase, nRef, nTest)
%randomly select k people from haplotpye sequence hapSeq, make sure that
%the 2 happlotype sequences of one individual is adjacent row
[m, nSnps] = size(hapSeq);
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));   %reset the random number generator
totalSample = nCase + nRef + nTest;
%if total number of sequence exceed the total number of indivuals
if 2*totalSample > m
    e = MException('signPower:notEnoughUniqueSeq', 'xzhou:not enough unique sequence');
    throw(e);
end
idx = randsample(floor(m/2), totalSample);
idxCase = sort([idx(1:nCase)*2; idx(1:nCase)*2-1]);
idxRef = sort([idx(nCase+1:nCase+nRef)*2; idx(nCase+1:nCase+nRef)*2-1]);
idxTest = sort([idx(nCase+nRef+1:end)*2; idx(nCase+nRef+1:end)*2-1]);

caseSeq = hapSeq(idxCase, :);
refSeq = hapSeq(idxRef, :);
testSeq = hapSeq(idxTest, :);
end