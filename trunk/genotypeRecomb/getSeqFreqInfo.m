function [freqInfo] = getSeqFreqInfo(seq, blocks)
% returns the frequency informatin of each block
[m ~] = size(blocks);
freqInfo = cell(m,1);
for i = 1:m
    freqInfo{i,1} = gBlockFreq(seq, blocks(i,:));
end
end