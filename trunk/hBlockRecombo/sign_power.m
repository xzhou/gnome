function [caseTr, refTr, testTr, caseTp, refTp, testTp, caseTrRefSign, testTrRefSign] = sign_power(hapSeq, nCase, nRef, nTest, trial, useEstR)
%explore the relationship between sign and power. Previously we though some
%sign has larger power. And even 10% sign has significant power. I think
%it's now stable since some pair has negative effects to the power,
%removing those sign or location might increase the sign. But the attacker
%can not tell the difference. 

%@hapSeq is the 0-1 encoding of hap snps sequence
%@caseSize is the population size of case and reference

%uniqueHapSeq = unique(hapSeq, 'rows');

[m, nSnps] = size(hapSeq);

%total number of individuals, each individual has 2 sequence
totalSample = nCase + nRef + nTest;

%if total number of sequence exceed the total number of indivuals
if 2*totalSample > m
    e = MException('signPower:notEnoughUniqueSeq', 'xzhou:not enough unique sequence');
    throw(e);
end

%idx of individuals
idx = randsample(m/2, totalSample);
idxCase = sort([idx(1:nCase)*2; idx(1:nCase)*2-1]);
idxRef = sort([idx(nCase+1:nCase+nRef)*2; idx(nCase+1:nCase+nRef)*2-1]);
idxTest = sort([idx(nCase+nRef+1:end)*2; idx(nCase+nRef+1:end)*2-1]);

caseSeq = hapSeq(idxCase, :);
refSeq = hapSeq(idxRef, :);
testSeq = hapSeq(idxTest, :);

%convert hapseq to genoseq
caseGenoSeq = combineHapSeq(caseSeq);
refGenoSeq = combineHapSeq(refSeq);
testGenoSeq = combineHapSeq(testSeq);

if useEstR == 0
    caseR = corrcoef(caseSeq);
    refR = corrcoef(refSeq);
else
    caseR = estimateR(caseGenoSeq);
    refR = estimateR(refGenoSeq);
end

caseR(isnan(caseR)) = 0;%invariant sites
refR(isnan(refR)) = 0;%invariant sites

caseP = sum(caseSeq)/nCase/2;
refP = sum(refSeq)/nRef/2;

nSeg = 10; %divide 1-100 into 10 segment
caseTr = zeros(nSeg, trial, 2*nCase);
refTr = zeros(nSeg, trial, 2*nRef);
testTr = zeros(nSeg, trial, 2*nTest);

caseTp = zeros(nCase, 1);
refTp = zeros(nRef, 1);
testTp = zeros(nTest, 1);

%calculate Homer's attack
%case
for k = 1:nCase
    Yk = caseGenoSeq(k, :);
    Tp = getTp(Yk, caseP, refP);
    caseTp(k) = Tp;
end
%ref
for k = 1:nRef
    Yk = refGenoSeq(k, :);
    Tp = getTp(Yk, caseP, refP);
    refTp(k) = Tp;
end
%test
for k = 1:nTest
    Yk = testGenoSeq(k, :);
    Tp = getTp(Yk, caseP, refP);
    testTp(k) = Tp;
end

fprintf(1, 'Homer: avg case = %f, ref = %f, test = %f\n', mean(caseTp), mean(refTp), mean(testTp));

%copy sign of reference r
caseTrRefSign = zeros(2*nCase, 1);
caseRRefSign = abs(caseR).*sign(refR);
for k = 1:2*nCase
    Yc = caseSeq(k, :);
    Tr = getTr(Yc, caseRRefSign, refR);
    caseTrRefSign(k) = Tr;
end
testTrRefSign = zeros(2*nTest, 1);
for k = 1:2*nTest
    Yc = testSeq(k, :);
    Tr = getTr(Yc, caseRRefSign, refR);
    testTrRefSign(k) = Tr;
end

%calculate power
%for different recover rate
for i = 1:nSeg
    %fprintf(1, '%d\n', i);
    p = 0.1 * i;
    fprintf(1, 'p = %f%%\n', p*100);
    %try multiple times
    for j = 1:trial
        pMask = getMaskP(nSnps, p);
        pCaseR = caseR .* pMask;
        %for each individual
        for k = 1:2*nCase
            Yc = caseSeq(k,:);
            Tr = getTr(Yc, pCaseR, refR);
            caseTr(i, j, k) = Tr;
        end
        for k = 1:2*nRef
            Yr = refSeq(k, :);
            Tr = getTr(Yr, pCaseR, refR);
            refTr(i, j, k) = Tr;
        end
        for k = 1:2*nTest
            Yt = testSeq(k, :);
            Tt = getTr(Yt, pCaseR, refR);
            testTr(i, j, k) = Tt;
        end
    end
end
end