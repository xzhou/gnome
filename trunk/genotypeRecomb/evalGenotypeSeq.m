function [currentQuality currentR blockMask] = evalGenotypeSeq(targetRs, currentGenotypeSeq, block1, block2, config)
    
    alpha = config.alpha;
    smallFilter = config.smallFilter;
    
    blockMask = getBlockMaskForEval(targetRs, block1, block2);
        
    newR = estimateR(currentGenotypeSeq);
    newRs = newR.*newR;
    targetRs = targetRs.*blockMask;
    newRs = newRs .* blockMask;
    
    %remove diagnal elements
    targetRs(logical(eye(size(targetRs)))) = 0;
    newRs(logical(eye(size(newRs)))) = 0;
    
    %filter the small r values
    targetRs(targetRs < smallFilter) = 0;
    newRs(newRs < smallFilter) = 0;
    
    currentQuality = nansum(nansum(abs(targetRs - newRs)))/2;
    currentR = newR;
    
end