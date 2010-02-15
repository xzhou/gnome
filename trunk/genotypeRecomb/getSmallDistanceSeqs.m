function [xyVal] = getSmallDistanceSeqs(refHapPool, nSeq, targetRs, alleleMapping)
%calculate the r square of a single sequence and then compare the distance
%to target rsquare, and put the smallest nSeq sequence to bestSeq, this is
%very similiar to homer's test

[nHap, nSnps] = size(refHapPool);

nGenoCombo = nHap*(nHap-1)/2;

rsdiff = NaN(nHap,nHap);

%for each possible haplotype combination, calculate the r square distance
for i = 1:nHap
    for j = i:nHap
        A = refHapPool(i,:);
        B = refHapPool(j,:);
        AB = haplotype2genotype([A;B], alleleMapping);
        rs = calcSelfCorr(AB);
        diff = getRsDiff(rs, targetRS);
        rsdiff(i,j) = diff;
        %rsdiff(j,i) = diff;
    end
end

if nGenoCombo > nSeq
    xy = topk(rsdiff, nSeq);
else
    %not enough combination, select all and padding according to frequency
    nMissing = nSeq - nGenoCombo;
    xyVal = topk(rsdiff, nGenoCombo);
    
    %get higher frequency conbination
    mxy = randomTopK(rsDiff, nMissing);
    xyVal = [xyVal; mxy];
end


end





