function [signAgreeRate] = blockSignRate(targetR, newR, blocks, newBlock)
    blockMask = getBlockMaskForEval(targetR, blocks, newBlock);
    
    targetR = targetR.*blockMask;
    newR = newR.*blockMask;
    
    signDiff = sign(targetR.*newR);
    signDiff(logical(eye(size(signDiff)))) = 0;
    
    
    nDiff = sum(sum(signDiff==-1))/2;
    totalSign = sum(sum(abs(signDiff)))/2;
    
    signAgreeRate = nDiff/totalSign;
    
end