function [value] = SignRate(targetR, sampleR)
    targetSign = sign(targetR);
    sampleSign = sign(sampleR);
    
    signDiff = targetSign.*sampleSign;
    signDiff(logical(eye(size(signDiff)))) = 0;
    
    correctSign = sum(sum(double(signDiff == 1)))/2;
    
    [m n] = size(targetR);
    
    totalSign = m*(m-1)/2;
    
    value = correctSign/totalSign;
end