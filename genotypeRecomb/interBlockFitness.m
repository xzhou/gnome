function rDifference = interBlockFitness(genotypeBlock, caseBlockSeq, block1, block2)

%[r1 PA1 count1] = estimateR(genotypeBlock);
r1 = calcR(genotypeBlock);

caseBlockSeq = haplotype2genotype(caseBlockSeq);
[r2 PA2 count2] = estimateR(caseBlockSeq);

blockMask = getBlockMaskForEval(caseBlockSeq, block1, block2);

r1s = r1.*r1;
r2s = r2.*r2;

r1s = r1s.*blockMask;
r2s = r2s.*blockMask;

r1s(logical(eye(size(r1s)))) = 0;
r2s(logical(eye(size(r2s)))) = 0;

rDifference = sum(sum(abs(r1s-r2s)))/2;