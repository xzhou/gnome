%evaluate normalized r square difference and p difference 
function [newQuality] = getQuality(targetRs, currentRs, targetFreq, currentFreq, alpha)
    %we think R square has more weight than single allele frequence
    if nargin == 4
        alpha = 0.6;
    end

    pDiff = abs(targetFreq - currentFreq);
    normalDiff = sum(pDiff)/length(targetFreq);   %normalize
    
    rDiff = calcRDiff(targetRs, currentRs);
    [m, n] = size(targetRs);
    if m ~= n
        e = MException('blockReconstruct:x', 'in consistent dimenstion');
        throw(e);
    end
    nElements = m*(m-1)/2;
    normalRDiff = rDiff/nElements;
    
    newQuality = alpha*normalRDiff + (1-alpha)*normalDiff;
    newQuality = 100*newQuality;
end

function [newQuality] = calcRDiff(targetRs, newRs, smallRFilter)
    
    if nargin == 2
        smallRFilter = 0.0;
    end
    %make sure remove the diagnal elements
    targetRs(logical(eye(size(targetRs)))) = 0;
    newRs(logical(eye(size(newRs)))) = 0;
   
    %filter those points whoes Rs is 0 in either Case or Ref
    targetRs(logical(targetRs==0))=NaN;
    newRs(logical(newRs==0))=NaN;
    
    %filter the small r values
    targetRs(targetRs < smallRFilter) = 0;
    newRs(targetRs < smallRFilter) = 0;
    
    newQuality = nansum(nansum(abs(targetRs - newRs)))/2;
end