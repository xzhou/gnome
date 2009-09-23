function [finalSeq finalR finalSignRate finalQual] = newHBRecombo(targetR, targetSeq, currentSeq, blocks, newBlock, alleleMapping, f)

    
    debugMode = 1;

    targetRs = targetR.*targetR;
    
    %TODO calculate the size
    expT = numel(targetR)/2*0.000001;
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
    fprintf(1, 'initQ = %f \tfinalQ = %f initSR = %f \tfinalSR = %f\n', initQuality, finalQual, initSignRate, finalSignRate);
    if nargin == 7
        fprintf(f, 'initQ = %f \tfinalQ = %f initSR = %f \tfinalSR = %f\n', initQuality, finalQual, initSignRate, finalSignRate);
    end
end