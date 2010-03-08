function [caseTr, refTr] = sign_power(hapSeq, nCase, nRef, trial)
%explore the relationship between sign and power. Previously we though some
%sign has larger power. And even 10% sign has significant power. I think
%it's now stable since some pair has negative effects to the power,
%removing those sign or location might increase the sign. But the attacker
%can not tell the difference. 

%@hapSeq is the 0-1 encoding of hap snps sequence
%@caseSize is the population size of case and reference

uniqueHapSeq = unique(hapSeq, 'rows');

[m, nSnps] = size(uniqueHapSeq);

totalSample = nCase + nRef;
if totalSample > m
    e = MException('signPower:notEnoughUniqueSeq', 'xzhou: not enough unique sequence');
    throw(e);
end

idx = randsample(m, totalSample);

caseSeq = uniqueHapSeq(idx(1:nCase), :);
refSeq = uniqueHapSeq(idx(nCase+1:end), :);

caseR = corrcoef(caseSeq);
refR = corrcoef(refSeq);

caseR(isnan(caseR)) = 0;%invariant sites
refR(isnan(refR)) = 0;%invariant sites

nSeg = 10; %divide 1-100 into 10 segment
caseTr = zeros(nSeg, trial, nCase);
refTr = zeros(nSeg, trial, nRef);

%calculate power
for i = 1:nSeg
    fprintf(1, '%d\n', i);
    p = 0.1 * i;
    for j = 1:trial
        pMask = getMaskP(nSnps, p);
        pCaseR = caseR .* pMask;
        for k = 1:nCase
            Yc = caseSeq(k,:);
            Tr = getTr(Yc, pCaseR, refR);
            caseTr(i, j, k) = Tr;
        end
        for k = 1:nRef
            Yr = refSeq(k, :);
            Tr = getTr(Yr, pCaseR, refR);
            refTr(i, j, k) = Tr;
        end
    end
end
end