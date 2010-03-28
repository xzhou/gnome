function rDifference = innerBlockFitness(genotypeBlock, tempCaseRs, tempCaseAlleleFreq, alleleMapping)

%[r1 PA1 count1] = estimateR(genotypeBlock);
currentR = calcR(genotypeBlock, alleleMapping);

% caseBlockSeq = haplotype2genotype(caseBlockSeq);
% [r2 PA2 count2] = estimateR(caseBlockSeq);
% r2 = calcR(caseBlockSeq);

currentRs = currentR.*currentR;

% r1s = r1.*r1;
% r2s = r2.*r2;

tempCaseR = tempCaseRs;

currentRs(logical(currentRs==0))=NaN;
tempCaseR(logical(tempCaseR==0))=NaN;

rDifference = nansum(nansum(abs(currentRs-tempCaseR)))/2;

currentAlleleFreq = GnomeCalculator.getSingleAlleleFreq(genotypeBlock, alleleMapping);

pDiff = abs(currentAlleleFreq - tempCaseAlleleFreq);
normalDiff = sum(pDiff)/length(currentAlleleFreq);   %normalize

alpha = 0.6;
[m, n] = size(tempCaseRs);
nElements = m*(m-1)/2;
normalRDiff = rDifference/nElements;

newQuality = alpha*normalRDiff + (1-alpha)*normalDiff;
rDifference = 100*newQuality;
