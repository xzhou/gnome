function [caseRef, coverRate, seq1Info, seq2Info] = analysisPhasedCaseRef(intSeq1, intSeq2, blocks)
%get hyplotype block cover information
[nseq, ~] = size(intSeq1);
seq1Info = getSeqFreqInfo(intSeq1, blocks);
seq2Info = getSeqFreqInfo(intSeq2, blocks);
[caseRef] = compareCaseRef(seq1Info, seq2Info, blocks);

[nBlocks, ~] = size(blocks);
coverRate.typeCoverate = 0.0;   %number of hyplotype shared
coverRate.seqCoverate = 0.0;    %number sequence shared
coverRate(nBlocks).typeCoverate = 0.0;
for i = 1:nBlocks
    [commonBlock, ~] = size(caseRef{i,1});
    [caseNBlock, ~] = size(seq1Info{i,1}.uniqueBlock);
    coverRate(i).typeCoverate = commonBlock*1.0/caseNBlock;
    if ~isempty(caseRef{i,1})
        coverRate(i).seqCoverate = sum(caseRef{i,1}(:,3))/nseq;
    else
        coverRate(i).seqCoverate = 0;
    end
end
end

function [freqInfo] = getSeqFreqInfo(seq, blocks)
% returns the frequency informatin of each block
[m ~] = size(blocks);
freqInfo = cell(m,1);
for i = 1:m
    freqInfo{i,1} = gBlockFreq(seq, blocks(i,:));
end
end