function [bestGenotypeCoverSeq] = selectBestCoverRef(refHapPool, caseSeq, wtccc1Conf)
%from the refHapPool, try to find the genotype that exists in caseSeq
%assuming that a sequence in case seq will reduce the r square difference.
%We will use a learning algorithm to get the best coverage

%configuration
[nHType, tmp] = size(refHapPool)
maxItr = nHType*(nHType-1)/2;

[targetRs, targetR, targetCounts] = estimateR(caseSeq);
[nSeq, nSnps] = size(caseSeq);

[initRefGSeq] = getSmallDistanceSeqs(refHapPool, nSeq);

itr = 0;
while itr < maxItr
    %find a new seq
    gSeq = drawOneGenotypeSeq(refHapPool);
    while containsSeq(caseSeq, seq)
        gSeq = drawOneGenotypeSeq(refHapPool);
    end
    %find the best replaceable seq in the caseSeq
    %scoreVector = [rowID, score, tested?];
    scoreVector = zeros(nSeq, 3);
    for i = nSeq
        caseSeqCopy = caseSeq;
        selectRow = aCaseSeq(i);
        testedRow = isTestedSeq(caseSeqCopy, selectRow, scoreVector(:,3));
        if testedRow <= 0
            caseSeqCopy(i,:) = gSeq;
            rs = estimateR(caseSeqCopy);
            rspDiff = calcRsqDiff(rs, p, targetRs, targetP);%find rs and p diff
            scoreVector(i) = [i, rspDiff, 1];
        else
            scoreVector(i) = [i, scoreVector(testedRow, 2), 1];
        end
    end
    %find the highest score
    [minScore, minIdx] = min(scoreVector(:,2));
    bestReplaceIdx = scoreVecotor(minIdx,1);
    caseSeq(bestReplaceIdx, :) = gSeq;
    scoreDiff = minScore;
end
end

function [existed] = containsSeq(seqs, seq)
existed = 0;
[nSeq, nSnps] = size(seqs);
for i = 1:nSeq
    if sum(seqs(i) ~= seq) == 0
        existed = i;
        break;
    end
end
end

function [tested] = isTestedSeq(seqs, seq, testBits)
len = length(testBits);
tested = 0;
for i = 1:len
    if testBits(i) == 1 && sum(seqs(i,:) ~= seq) == 1
        tested = i;
        break;
    end
end
end