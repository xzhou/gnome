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
        blocks(i, 4) = 1;
    end
    
    finalSeq = currentSeq;
    finalR = calcR(finalSeq);
end

function [nLearnedBlocks totalBlocks] = getLearnedBlocks(blocks)
    nLearnedBlocks = sum(blocks(:,4));
    [totalBlocks ncol] = size(blocks);
end

function [finalSeq finalR finalSignRate] = newHBRecombo(targetR, targetSeq, currentSeq, blocks, newBlock, alleleMapping)
    debugMode = 1;

    targetRs = targetR.*targetR;
    
    %TODO calculate the size
    expT = numel(targetR)/2*0.0001;
    maxIT = 10000;
    
    [currentQuality currentR] = blockEvaluateSeq(targetRs, currentSeq, alleleMapping, blocks, newBlock);
    initSignRate = blockSignRate(targetR, currentR, blocks, newBlock);
    
    initQuality = currentQuality;
    
    %TODO new SA algorithm
    t = 0;
    signRate = initSignRate;
    while t < maxIT
        t = t + 1;
        newSampleSeq = newMutate(currentSeq, newBlock);
        [newQuality newR] = blockEvaluateSeq(targetRs, newSampleSeq, alleleMapping, blocks, newBlock);
        Qdiff = newQuality - currentQuality;
        if Qdiff < 0
            p = 1;
        else
            p = 1/(1+exp(Qdiff/expT));
        end
        x = rand(1);
        if x < p
            currentSeq = newSampleSeq;
            currentQuality = newQuality;
            currentR = newR;
            signRate = blockSignRate(targetR, newR, blocks, newBlock);
            
            if debugMode == 1
                %fprintf(1, '%d currentQ = %f, signRate = %f\n', t, currentQuality, signRate);
            end
        else
            %debug information
        end
    end
    
    finalSeq = currentSeq;
    finalQual = currentQuality;
    finalR = currentR;
    finalSignRate = signRate;
    fprintf(1, 'initQ = %f \tfinalQ = %f initSR = %f \tfinalSR = %f\n', initQuality, finalQual, initSignRate, finalSignRate);
end

function [currentSeq] = newMutate(currentSeq, newBlock)
    [m n] = size(currentSeq);
    a = newBlock(1,1);  %start snps
    b = newBlock(1,2);  %end snps
    
    idx1 = randi([1 m]);
    idx2 = randi([1 m]);
    
    %mutate
    while currentSeq(idx1, a:b) == currentSeq(idx2, a:b)
        idx1 = randi([1 m]);
        idx2 = randi([1 m]);
    end
    
    %swap the sequence
    temp = currentSeq(idx1, a:b);
    currentSeq(idx1, a:b) = currentSeq(idx2, a:b);
    currentSeq(idx2, a:b) = temp;
end

