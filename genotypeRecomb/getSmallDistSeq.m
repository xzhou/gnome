function [gSeq] = getSmallDistSeq(refHapPool, nSeq, targetRs, targetF, alleleMapping, blocks)
%this function try to find a optimal reference sequences that best cover
%the target sequence. The heuristic here is that if a genotype sequence is
%in the target sequence, it's distance to the targetRs and target F should
%be statically smaller than those that is not in target sequence. 
nBlock = length(refHapPool);
gSeq = [];

for i = 1:nBlock
    a = blocks(i, 1);
    b = blocks(i, 2);
    blockRefPool = refHapPool(i).uniqueBlock;
    blockRs = targetRs(a:b, a:b);
    blockF = targetF(a:b);
    blockSeq = getSmallDistanceBlock(blockRefPool, nSeq, blockRs, blockF, alleleMapping(a:b));
    gSeq = [gSeq, blockSeq];
end

end