function [caseTr, refTr] = sign_power(hapSeq, nCase, nRef, trial)
%explore the relationship between sign and power. Previously we though some
%sign has larger power. And even 10% sign has significant power. I think
%it's now stable since some pair has negative effects to the power,
%removing those sign or location might increase the sign. But the attacker
%can not tell the difference. 

%@hapSeq is the 0-1 encoding of hap snps sequence
%@caseSize is the population size of case and reference

uniqueHapSeq = unique(hapSeq, 'rows');

[m, n] = size(uniqueHapSeq);

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
caseTr = zeros(nSeg, nCase, trial);
refTr = zeros(nSeg, nRef, trial);

%calculate power
for i = 1:nSeg
    fprintf(1, '%d\n', i);
    p = 0.1 * i;
    for j = 1:nCase
        Y = caseSeq(j,:);
        caseTr(i,j, :) = getTr_P_N(Y, caseR, refR, p, trial);
    end
    for j = 1:nRef
        Y = refSeq(j,:);
        refTr(i,j, :) = getTr_P_N(Y, caseR, refR, p, trial);
    end
end

end

function [Tr] = getTr_P_N(Y, caseR, refR, p, trial)
%calculate the power of Y using 100p percent signs, try trial number of
%times
Tr = zeros(trial, 1);
for i = 1:trial
    Tr(i) = getTr_P(Y, caseR, refR, p);
end
end

function [tr] = getTr_P(Y, caseR, refR, p)
%randomly inverse 1-p signs
[m n] = size(caseR);
mask = ones(m, n);
for i = 1:m-1
    x = (rand(m - i, 1) < p)*2 - 1;
    mask(i, i+1:m) = x;
end
mask = copyUpperToLower(mask);
caseR = caseR .* mask;
tr = getTr(Y, caseR, refR);
tr = tr/sqrt(n*(n-1)/2);
end