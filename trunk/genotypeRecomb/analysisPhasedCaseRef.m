function [caseRef, coverRate, seq1Info, seq2Info] = analysisPhasedCaseRef(intSeq1, intSeq2, blocks)
    [nseq, tmp] = size(intSeq1);
    seq1Info = getSeqFreqInfo(intSeq1, blocks);
    seq2Info = getSeqFreqInfo(intSeq2, blocks);
    [caseRef] = compareCaseRef(seq1Info, seq2Info, blocks);
    
    [nBlocks, n] = size(blocks);
    coverRate.typeCoverate = 0.0;   %number of hyplotype shared
    coverRate.seqCoverate = 0.0;    %number sequence shared
    coverRate(3).typeCoverate = 0.0;
    for i = 1:nBlocks
        [commonBlock, tmp] = size(caseRef{i,1});
        [caseNBlock, tmp] = size(seq1Info{i,1}.uniqueBlock);
        coverRate(i).typeCoverate = commonBlock*1.0/caseNBlock;
        coverRate(i).seqCoverate = sum(caseRef{1,1}(:,3))/nseq;
    end
end

function [freqInfo] = getSeqFreqInfo(seq, blocks)
% returns the frequency informatin of each block
    [m n] = size(blocks);
    freqInfo = cell(m,1);
    for i = 1:m
        freqInfo{i,1} = gBlockFreq(seq, blocks(i,:));
    end
end