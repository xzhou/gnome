function [sampledSeq] = blockSampleHapSeq(hapSeqNoID, config)
%%BLOCK_SAMPLE_HAP_SEQ will sample each block independently and assemble
%%them into a haplotype seq
genoSeqSize = config.caseSize;
[nSeqPool, nSnpPool] = size(hapSeqNoID);

blocks = config.blocks;
[nBlock, tmp] = size(blocks);
sampledSeq = [];
for i = 1:nBlock
    a = blocks(i,1);
    b = blocks(i,2);
    oneBlock = hapSeqNoID(:,a:b);
    newBlock = mySample(oneBlock, config);
    sampledSeq = [sampledSeq, newBlock];
end

end

function [sampledHapSeq] = mySample(hapSeqPool, config)
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