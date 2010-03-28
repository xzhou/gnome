function rDifference = interBlockFitness(genotypeBlock, caseBlockSeq, block1, block2)

%[r1 PA1 count1] = estimateR(genotypeBlock);
r1 = calcR(genotypeBlock);

% caseBlockSeq = haplotype2genotype(caseBlockSeq);
% [r2 PA2 count2] = estimateR(caseBlockSeq);

r2 = calcR(caseBlockSeq);

blockMask = getBlockMaskForEval(r2, block1, block2);

r1s = r1.*r1;
r2s = r2.*r2;

r1s = r1s.*blockMask;
r2s = r2s.*blockMask;

%filter those points whoes Rs is 0 in either Case or Ref
r1s(logical(r1s==0))=NaN;
r2s(logical(r2s==0))=NaN;

%make sure remove the diagnal elements
r1s(logical(eye(size(r1s)))) = 0;
r2s(logical(eye(size(r2s)))) = 0;

rDifference = nansum(nansum(abs(r1s-r2s)))/2;