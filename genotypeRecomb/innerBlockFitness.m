function rDifference = innerBlockFitness(genotypeBlock, caseBlockSeq)

%[r1 PA1 count1] = estimateR(genotypeBlock);
r1 = calcR(genotypeBlock);

% caseBlockSeq = haplotype2genotype(caseBlockSeq);
% [r2 PA2 count2] = estimateR(caseBlockSeq);
r2 = calcR(caseBlockSeq);

r1s = r1.*r1;
r2s = r2.*r2;

rDifference = sum(sum(abs(r1s-r2s)))/2;
