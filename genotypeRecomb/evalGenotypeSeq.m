function [currentQuality currentR] = evalGenotypeSeq(targetRs, currentGenotypeSeq, alleleMapping, block1, block2, alleleMapping, alpha, smallFilter)
    if nargin <= 5
        alpha = 0.1;
        smallFilter = 0;
    end
    
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
    newRs(newRs < samllFilter) = 0;
    
    newQuality = sum(sum(abs(targetRs - newRs)))/2;
    
end