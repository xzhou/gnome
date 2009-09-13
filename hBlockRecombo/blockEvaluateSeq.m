function [newQuality newR] = blockEvaluateSeq(targetRs, newSampleSeq, alleleMapping, blocks, newBlock)
    %this is a new evaluation function that only count the
    %learnedSnpsBlocks and newBlock r value
    blockMask = getBlockMaskForEval(targetRs, blocks, newBlock);
    
    newR = calcR(newSampleSeq, alleleMapping);
    newRs = newR.*newR;
    
    %filter to the blocks
    targetRs = targetRs.*blockMask;
    newRs = newRs.*blockMask;
    
    %make sure remove the diagnal elements
    targetRs(logical(eye(size(targetRs)))) = 0;
    newRs(logical(eye(size(newRs)))) = 0;
    
    newQuality = sum(sum(abs(targetRs - newRs)))/2;
end

