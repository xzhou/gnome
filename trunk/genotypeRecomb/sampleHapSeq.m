function [sampledHapSeq] = sampleHapSeq(hapSeqPool, config)
%hapSeqPool is a pool of hyplotype sequence
    genoSeqSize = config.caseSize;
    [nSeqPool, nSnpPool] = size(hapSeqPool);
    
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
    
    %sampling
    sampledIdx = zeros(2*genoSeqSize, 1);
    for i = 1:genoSeqSize
        idx1 = randi(nSeqPool);
        idx2 = randi(nSeqPool);
        sampledIdx([2*i-1,2*i]) = [idx1, idx2];
    end
    sampledHapSeq = hapSeqPool(sampledIdx, :);
end