function [finalR, finalRDiff, finalQuality] = recursiveLearning(caseSeq4, refSeq4, blocks)
    
    debugMode = 1;
    %calculation
    alleleMapping = getMajorAllele(refSeq4);
    refR = calcR(refSeq4, alleleMapping);
    caseR = calcR(caseSeq4, alleleMapping);
    
    %add the size of the block to the blocks
    [nBlock nBlockCol] = size(blocks);
    blocks = [blocks blocks(:,2)-blocks(:,1) + 1 zeros(nBlock, 1)];
    
    firstBlock = getNextBlock(blocks);
    blocks(firstBlock, 4) = 1;
    
    %learn from current Sequence
    currentSeq = refSeq4;
    
    while any(blocks(:,4)~=1)
        %get next block row index
        i = getNextBlock(blocks);
        
        newStart = blocks(i, 1);
        newEnd = blocks(i, 2);
        
        newBlock = blocks(i,:);
        
        %TODO: inside block learning
        
        
        %TODO: we need a more sophiscated algorithm to detecth the
        %convergence of the learning
        finalSignRate = 0;
        
        nLearnedBlocks = getLearnedBlocks(blocks);
        fprintf(1, 'learning block %d, %d blocks has fixed\n', i, nLearnedBlocks);
        
        [tempFinalSeq tempFinalR tempFinalSignRate] = newHBRecombo(caseR, caseSeq4, currentSeq, blocks, newBlock, alleleMapping);
        currentSeq = tempFinalSeq;
        currentR = tempFinalR;
        diff = caseSeq4 - currentSeq;
        blocks(i, 4) = 1;
    end
    
    finalSeq = currentSeq;
    finalR = calcR(finalSeq);
end





