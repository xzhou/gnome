function [newQuality newR blockMask] = blockEvaluateSeq(targetRs, newSampleSeq, refSeq4, alleleMapping, blocks, newBlock, smallRFilter, alpha)
    %this is a new evaluation function that only count the
    %learnedSnpsBlocks and newBlock r value
    if nargin <= 5
        smallRFilter = 0;
    end
    
    if nargin <= 7
        alpha = 0.8;
    end
    
    blockMask = getBlockMaskForEval(targetRs, blocks, newBlock);
    
    newR = calcR(newSampleSeq, alleleMapping);
    newRs = newR.*newR;
    
    %filter to the blocks
    targetRs = targetRs.*blockMask;
    newRs = newRs.*blockMask;
    
    %filter those points whoes Rs is 0 in either Case or Ref
    targetRs(logical(targetRs==0))=NaN;
    newRs(logical(newRs==0))=NaN;
    
    %make sure remove the diagnal elements
    targetRs(logical(eye(size(targetRs)))) = 0;
    newRs(logical(eye(size(newRs)))) = 0;
    
    %filter the small r values
    targetRs(targetRs < smallRFilter) = 0;
    newRs(targetRs < smallRFilter) = 0;
    
    rDiff = nansum(nansum(abs(targetRs - newRs)))/2;
    
    [m, n] = size(targetRs);
    if m ~= n
        e = MException('blockReconstruct:x', 'in consistent dimenstion');
        throw(e);
    end
    nElements = blocks(1,3).*newBlock(1,3);
    normalRDiff = rDiff/nElements;
    
    %measure the C00 distance between ref and currentSeq
%     [newR newC00 newC01 newC10 newC11] = calcPairwiseFreq(newSampleSeq, alleleMapping);
%     [refR refC00 refC01 refC10 refC11] = calcPairwiseFreq(refSeq4, alleleMapping);
%     
%     c00Diff = sum(sum(abs(newC00.*blockMask-refC00.*blockMask)));
%     normalC00Diff = c00Diff/nElements/500;
%     
%     newQuality = normalRDiff*alpha + (1-alpha)*normalC00Diff; 
    newQuality = normalRDiff;
    
end

