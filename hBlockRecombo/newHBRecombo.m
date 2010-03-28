function [finalSeq finalR finalSignRate finalQual blockMask] = newHBRecombo(targetR, targetSeq, refSeq4, refC00, currentSeq, blocks, newBlock, alleleMapping, smallFilter)

    debugMode = true;

    targetRs = targetR.*targetR;
    
    %TODO calculate the size
    expT = numel(targetR)/2*1e-10;
    maxIT = 10000;
    
    %For fast computing
    blockMask = getBlockMaskForEval(targetRs, blocks, newBlock);
    targetRs = targetRs.*blockMask;
    
    currentR = calcR(currentSeq, alleleMapping);
    
    NaNmask = and((targetRs~=0), (currentR~=0));
    targetRs = targetRs.*NaNmask;
    
    AllMask = blockMask.*NaNmask;
    
    [currentQuality currentR] = blockEvaluateSeq(targetRs, currentSeq, refSeq4, refC00, alleleMapping, blocks, newBlock, AllMask, smallFilter);
    initSignRate = blockSignRate(targetR, currentR, blocks, newBlock);
    
    initQuality = currentQuality;
    
    %TODO new SA algorithm
    t = 0;
    signRate = initSignRate;
    
    while t < maxIT
        t = t + 1;
        newSampleSeq = newMutate(currentSeq, newBlock);
        [newQuality newR] = blockEvaluateSeq(targetRs, newSampleSeq, refSeq4, refC00, alleleMapping, blocks, newBlock, AllMask, smallFilter);
        Qdiff = newQuality - currentQuality;
        
        %Some problem here with this algorithm??
        if Qdiff < 0
            p = 1/(1+exp(Qdiff/expT));
            %p = 1;
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
    if debugMode
        fprintf(1, 'initQ = %.20f \tfinalQ = %.20f initSR = %f \tfinalSR = %f\n', initQuality, finalQual, initSignRate, finalSignRate);
    end
end