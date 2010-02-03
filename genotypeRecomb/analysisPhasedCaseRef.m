function [caseRef, coverRate] = analysisPhasedCaseRef(intSeq1, intSeq2, blocks)
    seq1Info = getSeqFreqInfo(intSeq1, blocks);
    seq2Info = getSeqFreqInfo(intSeq2, blocks);
    [caseRef] = compareCaseRef(seq1Info, seq2Info, blocks);
        
    [nBlocks, n] = size(blocks);
    coverRate = zeros(nBlocks, 1);
    for i = 1:nBlocks
        [commonBlock, tmp] = size(caseRef{i,1});
        [caseNBlock, tmp] = size(seq1Info{i,1}.uniqueBlock);
        coverRate(i,1) = commonBlock*1.0/caseNBlock;
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