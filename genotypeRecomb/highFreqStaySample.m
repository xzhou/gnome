function [sampledSeq] = highFreqStaySample(refPool, caseSize, T, blocks)
%%for each block, we keep the high frequency sequencies, and random sample
%%low frequencies blocks. T is the threshhold of frequency

[nRefSeq, ~] = size(refPool);
ratio = nRefSeq/caseSize;
[nBlock, ~] = size(blocks);

refFreqInfo = getSeqFreqInfo(refPool, blocks);

%keep high frequency
sampledSeq = [];
for i = 1:nBlock
    blocki = [];
    uniqueBlock = refFreqInfo{i}.uniqueBlock;
    freq = refFreqInfo{i}.freq;
    nUniqueBlock = length(freq);
    for j = 1:nUniqueBlock
        if freq(j) > T
            nkeep = round(freq(j)/ratio);
            blocki = [blocki;repmat(uniqueBlock(j), nkeep, 1)];
            %remove if keep
            uniqueBlock(j) = [];
            freq(j) = [];
        end
    end
    
    %sample missing sequencies from low freq seq
    [nHighFreq ~] = size(blocki);
    nMissing = caseSize - nHighFreq;
    missingSeq = distSampleK(uniqueBlock, freq, nMissing);
    blocki = [blocki, missingSeq];
    sampledSeq = [sampledSeq, blocki];
end
end