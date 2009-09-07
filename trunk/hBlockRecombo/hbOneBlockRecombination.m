function [finalSeq, finalR, finalSignRate] = hbOneBlockRecombination(targetR, refSeq, alleleMapping, blocks)
%hbOneBlockRecombination will use haplotype recombination to recover sign
    %config
    expT = numel(targetR)/2*0.0001;
    maxIT = 100000;
    
    %get the blocks
    a1 = blocks(1,1);
    b1 = blocks(1,2);
    a2 = blocks(2,1);
    b2 = blocks(2,2);
    
    block1Size = b1 - a1 + 1;
    block2Size = b2 - a2 + 1;
    
    refR = calcR(refSeq, alleleMapping);
    
    [h1 f1] = getHaplotype(refSeq, 1:block1Size);
    [h2 f2] = getHaplotype(refSeq, (block1Size+1):(block1Size+block2Size));
    
    currentSeq = refSeq;
    
    targetRs = targetR.*targetR;
    
    [currentQuality currentR] = evaluateSeqAbs(targetRs, currentSeq, alleleMapping);
    
    initSignRate = SignRate(targetR, currentR);
    initQuality = currentQuality;
    
    t = 0;
    signRate = initSignRate;
    while t < maxIT
        t = t+1;
        newSampleSeq = mutate(currentSeq, block1Size);
        [newQuality newR] = evaluateSeqAbs(targetRs, newSampleSeq, alleleMapping);
        Qdiff = newQuality - currentQuality;
        if Qdiff < 0
            %p = 1;
            p = 1/(1+exp(Qdiff/expT));
            %fprintf(1, 'BINGO! p = %f\n', p);
        else
            p = 1/(1+exp(Qdiff/expT));
            %p = 0;
        end
        x = rand(1);
        if x < p
            currentSeq = newSampleSeq;
            currentQuality = newQuality;
            currentR = newR;
            signRate = SignRate(targetR, newR, block1Size);
            %fprintf(1, 'Q = %f\t signRate = %f\n', currentQuality, signRate);
        else
            %null
        end
    end
    
    finalSeq = currentSeq;
    finalQual = currentQuality;
    finalR = currentR;
    finalSignRate = signRate;
    fprintf(1, 'initQ = %f \tfinalQ = %f initSR = %f \tfinalSR = %f\n', initQuality, finalQual, initSignRate, finalSignRate);
    SvsQ = [finalSignRate finalQual];
end