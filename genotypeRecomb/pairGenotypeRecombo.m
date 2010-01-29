function [finalSeq finalR fianlSignRate finalQual blockMask] = pairGenotypeRecombo(targetR, caseSeq, startingSeq, block1, block2, config)
    %function pairGenotypeRecombo recombine sequences
    %======= config =========
    expT = numel(targetR)/2*1e-10;
    maxIT = config.maxIT;
    verbose = config.verbose;
    %======= end config =====
    
    targetRs = targetR.*targetR;    
    %initialize
    [currentQuality currentR] = evalGenotypeSeq(targetRs, startingSeq, block1, block2, config);
    initSignRate = blockSignRate(targetR, currentR, block1, block2);
    initQuality = currentQuality;
    
    t = 0;
    signRate = initSignRate;
    currentSeq = startingSeq;
    while t < maxIT
        t = t + 1;
        newSeq = newMutate(currentSeq, block2);
        [newQuality newR blockMask] = evalGenotypeSeq(targetRs, newSeq, block1, block2, config);
        Qdiff = newQuality - currentQuality;
        if Qdiff < 0
            p = 1/(1+exp(Qdiff/expT));
        else
            p = 1/(1+exp(Qdiff/expT));
        end
        x = rand(1);
        if x < p
            currentSeq = newSeq;
            currentQuality = newQuality;
            currentR = newR;
            signRate = blockSignRate(targetR, newR, block1, block2);
            if verbose
                
            end
        else
           if verbose
                
           end
        end
    end
    
    finalSeq = currentSeq;
    finalQual = currentQuality;
    finalR = currentR;
    fianlSignRate = signRate;

end