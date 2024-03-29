function [result] = innerBlockLearning(caseBlock, caseFreq, refBlock, refFreq, refcaseFreq, alleleMapping, alpha)
    %alpha is the weight of r square in the learing process
    if nargin == 5
        alleleMapping = getMajorAllele(refSeq4);
        alpha = 0.6;
    elseif nargin == 6
        alpha = 0.6;
    end
    
    
    %defines the randomness of searching strategy
    expT = 1.0e-5;
    
    %% learning algorithm
    maxIt = 1e4;
    itr = 0;
    
    targetSeq = blockReconstruct(caseBlock, caseFreq);
    caseR = calcR(targetSeq, alleleMapping);
    caseRs = caseR.*caseR;
    caseAlleleFreq = GnomeCalculator.getSingleAlleleFreq(targetSeq, alleleMapping);

    currentSeq = blockReconstruct(refBlock, refFreq);
    refR = calcR(currentSeq, alleleMapping);
    currentR = refR;
    currentRs = currentR.*currentR;
    currentAlleleFreq = GnomeCalculator.getSingleAlleleFreq(currentSeq, alleleMapping);
    
    %initialize
    currentFreq = refFreq;
    currentQuality = getQuality(caseRs, currentR.*currentR, caseAlleleFreq, currentAlleleFreq, alpha);
    
    %fprintf(1, '%f\n', currentQuality);
    
    historySize = 1000;
    previousKQuality = zeros(historySize, 1);
    
    while itr < maxIt
        itr = itr + 1;
        [newSeq4 newFreq] = getNextSeq(refBlock, currentFreq);
        newR = calcR(newSeq4, alleleMapping);
        newRs = newR.*newR;
        %newP is the single allele frequency of the new allele frequency
        newP = GnomeCalculator.getSingleAlleleFreq(newSeq4, alleleMapping);
        newQuality = getQuality(caseRs, newRs, caseAlleleFreq, newP, alpha);
        Qdiff = newQuality - currentQuality;
        %using statistic hill climbing algorithm to change
        if Qdiff < 0
            p = 1/(1+exp(Qdiff/expT));
        else
            p = 1/(1+exp(Qdiff/expT));
        end
        x = rand(1);
        if x < p
            currentFreq = newFreq;
            currentQuality = newQuality;
            currentRs = newRs;
            currentR = newR;
            currentP = newP;
            %fprintf(1, 'itr = %d\t newQuality = %.15f\n', itr, newQuality);
        else
            %do nothing
        end
    end
    
    %end of learning
    result.finalFreq = currentFreq;
    result.finalQual = currentQuality;
    result.finalR = currentR;
    result.refFreq = refFreq;
    result.refcaseFreq = refcaseFreq;
    result.refBlock = refBlock;
    
    result.initSignRate = InnerBlockHelp.calcSignRate(caseR, refR);
    result.finalSignRate = InnerBlockHelp.calcSignRate(caseR, currentR);

    %calculate the frequency distance
    a = sum(abs(result.refcaseFreq - result.finalFreq));
    
    b = sum(abs(result.refcaseFreq - result.refFreq));
    
    %if fDistance < 1, the learning process is good
    result.fDistance = a*1.0/b;
end

function [newSeq, newBlockFreq] = getNextSeq(hyplotypes, currentBlockFreq)
    newBlockFreq = blockNaiveMutate(currentBlockFreq);
    newSeq = blockReconstruct(hyplotypes, newBlockFreq);
end


