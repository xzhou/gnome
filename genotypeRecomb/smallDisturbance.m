function [seq] = smallDisturbance(caseFreqInfo)
%keep different seqs types or rows in case sequence and resample the rest
%sequence
nBlock = length(caseFreqInfo);
seq = [];
for i = 1:nBlock
    blocki = [];
    uniqueSeq = caseFreqInfo{i,1}.uniqueBlock;
    freq = caseFreqInfo{i,1}.freq;
    nType = length(freq);
    nSeq = sum(freq);
    blocki = [uniqueSeq; distSampleK(uniqueSeq, freq, nSeq - nType)];
    seq = [seq, blocki];
end
end