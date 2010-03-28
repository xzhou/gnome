function [newQuality newR] = blockEvaluateSeq(targetRs, newSampleSeq, refSeq4, refC00, alleleMapping, blocks, newBlock, AllMask, smallRFilter, alpha)
    %this is a new evaluation function that only count the
    %learnedSnpsBlocks and newBlock r value
%     if nargin <= 5
%         smallRFilter = 0;
%     end
%     
    if nargin <= 9
        alpha = 0.8;
    end
    
    %blockMask = getBlockMaskForEval(targetRs, blocks, newBlock);
    
    if(blocks(1,1)<newBlock(1,2))
        block1 = blocks;
        block2 = newBlock;
    else
        block1 = newBlock;
        block2 = blocks;
    end
    
    [nS Len] = size(newSampleSeq);
    
    tempPartNewSeq(1:nS, 1:block1(1,3)) = newSampleSeq(1:nS, block1(1,1):block1(1,2));
    tempPartNewSeq(1:nS, (block1(1,3)+1):(block1(1,3)+block2(1,3))) = newSampleSeq(1:nS, block2(1,1):block2(1,2));

    tempAlleleMapping(:, 1:block1(1,3)) = alleleMapping(:, block1(1,1):block1(1,2));
    tempAlleleMapping(:, block1(1,3)+1:block1(1,3)+block2(1,3)) = alleleMapping(:, block2(1,1):block2(1,2));
    
    
    tempNewR = calcR(tempPartNewSeq, tempAlleleMapping);
    %tempNewRs = tempNewR.*tempNewR;
    newR = zeros(Len);
    newR(block2(1,1):block2(1,2), block1(1,1):block1(1,2)) = tempNewR((block1(1,3)+1):(block1(1,3)+block2(1,3)), 1:block1(1,3));
    newR(block1(1,1):block1(1,2), block2(1,1):block2(1,2)) = tempNewR(1:block1(1,3), (block1(1,3)+1):(block1(1,3)+block2(1,3)));
    newRs = newR.*newR;
    
    
%     newR = calcR(newSampleSeq, alleleMapping);
%     newRs = newR.*newR;
%     newRs = newRs.*blockMask;
    
    
    %filter to the blocks
    %targetRs = targetRs.*blockMask;
    newRs = newRs.*AllMask;
    
    %filter those points whoes Rs is 0 in either Case or Ref
    %NaNmask = (targetRs~=0).*(newRs~=0); 
    
    %NaNmask = and((targetRs~=0), (newRs~=0));
    %targetRs = targetRs.*NaNmask;
	%newRs = newRs.*NaNmask;
%     targetRs(logical(targetRs==0))=NaN;
%     newRs(logical(newRs==0))=NaN;
    
    %make sure remove the diagnal elements
    
    %targetRs(logical(eye(size(targetRs)))) = 0;
    %newRs(logical(eye(size(newRs)))) = 0;
    
    %filter the small r values
    %targetRs(targetRs < smallRFilter) = 0;
    %newRs(targetRs < smallRFilter) = 0;
    
    rDiff = sum(sum(abs(targetRs - newRs)))/2;
    
    [m, n] = size(targetRs);
    if m ~= n
        e = MException('blockReconstruct:x', 'in consistent dimenstion');
        throw(e);
    end
    nElements = blocks(1,3).*newBlock(1,3);
    normalRDiff = rDiff/nElements;
    
    %measure the C00 distance between ref and currentSeq
    
    newC00temp = calcC00(tempPartNewSeq, tempAlleleMapping);
    %     [newR newC001 newC01 newC10 newC11] = calcPairwiseFreq(newSampleSeq, alleleMapping);
    %     [refR refC00 refC01 refC10 refC11] = calcPairwiseFreq(refSeq4, alleleMapping);
    newC00 = zeros(Len);
    newC00(block2(1,1):block2(1,2), block1(1,1):block1(1,2)) = newC00temp((block1(1,3)+1):(block1(1,3)+block2(1,3)), 1:block1(1,3));
    newC00(block1(1,1):block1(1,2), block2(1,1):block2(1,2)) = newC00temp(1:block1(1,3), (block1(1,3)+1):(block1(1,3)+block2(1,3)));
    %
    c00Diff = sum(sum(abs(newC00.*AllMask-refC00.*AllMask)));
    normalC00Diff = c00Diff/nElements/500;
%     
     newQuality = normalRDiff*alpha + (1-alpha)*normalC00Diff; 
   % newQuality = normalRDiff;
    
end

