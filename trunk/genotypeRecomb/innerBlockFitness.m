function rDifference = innerBlockFitness(genotypeBlock, caseBlockSeq)

[r1 PA1 count1] = estimateR(genotypeBlock);

[r2 PA2 count2] = estimateR(caseBlockSeq);

rDifference = sum(sum(abs(r1-r2)));
