function [finalSeq finalR finalSignRate finalQual blockMask] = pairGenotypeRecombo(targetR, caseSeq, startingSeq, block1, block2, alpha)
    %function pairGenotypeRecombo recombine sequences
    %======= config =========
    expT = numel(targetR)/2*1e-10;
    maxIT = 1e4;
    smallFilter = 0.01;
    %======= end config =====
    
    verbose = true;
    targetRs = targetR.*targetR;
    
    %initialize
    [currentQuality currentR] = evalGenotypeSeq(targetRs, startingSeq, alleleMapping, block1, block2, alpha, smallFilter);
    initSignRate = blockSignRate(targetR, currentR, block1, block2);
    initQuality = currentQuality;
    
    itr = 0;
    signRate = initSignRate;
    currentSeq = startingSeq;
    while t < maxIT
        t = t + 1;
        newSampelSeq = newMutate(currentSeq, block2);
        [newQuality newR blockMask] = evalGenotypeSeq(targetRs, currentSeq, alleleMapping, block1, block2, alpha, smallFilter);
        Qdiff = newQuality - currentQuality;
        if Qdiff < 0
            p = 1/(1+exp(Qdiff/expT));
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
                %print
            end
        else
           if debugMode
               %print
           end
        end
    end
    
    finalSeq = currentSeq;
    finalQual = currentQuality;
    finalR = currentR;
    fianlSignRate = signRate;

end