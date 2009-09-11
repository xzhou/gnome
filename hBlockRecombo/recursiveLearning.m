function [finalR, finalRDiff, final] = recursiveLearning(caseSeq4, refSeq4, blockStruct)
    %calculation
    alleleMapping = getMajorAllele(refSeq4);
    refR = calcR(refSeq, alleleMapping);
    caseR = calcR(caseSeql, alleleMapping);
    
    %add the size of the block to the blocks
    blocks = [blocks blocks(:,1)-blocks(:,2), 0];
    
    %add the current status of the blocks in the 4th column, 0 means
    %unlearned
    blocks = [blocks zeros(length(blocks(:,1)), 1)];
    
    %continue if we have more unlearned blocks
    
    %learn from the beginning 
    currentSeq = refSeq4;
    
    while any(blocks(:,4)~=1)
        nextBlock = getNextBlock(blocks);
        
        learnedSnpsIndex = getAllIndex(blocks);
        newStart = blocks(i, 1);
        newEnd = blocks(i, 2);
        
        newIndex = newStart:newEnd;
        
        
        
        blocks(nextBlock, 4) = 1;
        
        

    end
end

