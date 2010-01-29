function [signAgreeRate] = blockSignRate(targetR, newR, blocks, newBlock)
    blockMask = getBlockMaskForEval(targetR, blocks, newBlock);
    
    targetR = targetR.*blockMask;
    newR = newR.*blockMask;
    
    signDiff = sign(targetR.*newR);
    signDiff(logical(eye(size(signDiff)))) = 0;
    
    
    nAgree = nansum(nansum(signDiff==1))/2;
    totalSign = nansum(nansum(abs(signDiff)))/2;
    
    signAgreeRate = nAgree/totalSign;
end