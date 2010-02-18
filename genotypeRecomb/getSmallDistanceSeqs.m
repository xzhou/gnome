function [gSeq, xyVal] = getSmallDistanceSeqs(refHapPool, nSeq, targetRs, targetF, alleleMapping)
%calculate the r square of a single sequence and then compare the distance
%to target rsquare, and put the smallest nSeq sequence to bestSeq, this is
%very similiar to homer's test

[nHap, nSnps] = size(refHapPool);

nGenoCombo = nHap*(nHap-1)/2;

rsdiff = NaN(nHap,nHap);

%for each possible haplotype combination, calculate the r square distance
for i = 1:nHap
    A = refHapPool(i,:);
    parfor j = i:nHap
        B = refHapPool(j,:);
        rs = calcSelfRs([A;B], alleleMapping);
        f = singleAlleleFreqOfHapSeq([A;B], alleleMapping);
        %AB = haplotype2genotype([A;B], alleleMapping);
        diff = getRsFDiff(rs, f, targetRs, targetF);
        rsdiff(i,j) = diff;
        %rsdiff(j,i) = diff;
    end
end

rsdiff = copyUpperToLower(rsdiff);

if nGenoCombo > nSeq
    xyVal = topk(rsdiff, nSeq);
else
    %not enough combination, select all and padding according to frequency
    nMissing = nSeq - nGenoCombo;
    xyVal = topk(rsdiff, nGenoCombo);
    
    %get higher frequency conbination
    mxy = randomTopK(rsDiff, nMissing);
    xyVal = [xyVal; mxy];
end

[m n] = size(xyVal);
if m ~= nSeq
    e = MException('getSmallDistanceSeqs:nSeq', 'please check xyVal');
    throw(e);
end

%construct the genotype
gSeq = constructGenotypeSeq(refHapPool, xyVal, alleleMapping);
end

function [gSeq] = constructGenotypeSeq(hapPool, indexes, alleleMapping)
%for all the indexes in indexes, select the haplotype sequence from hapPool
%and construct the genotype
[m n] = size(indexes);
[nHap, nSnps] = size(hapPool);
gSeq = zeros(m, nSnps);
for i = 1:m
    x = indexes(i, 1);
    y = indexes(i, 2);
    hseq = hapPool([x, y], :);
    oneGenoSeq = haplotype2genotype(hseq, alleleMapping);
    gSeq(i,:) = oneGenoSeq;
end
end





