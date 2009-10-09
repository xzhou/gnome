function [result] = innerBlockLearning(caseBlock, caseFreq, refBlock, refFreq, refcaseFreq, alleleMapping, alpha)
    if nargin == 5
        alleleMapping = getMajorAllele(refSeq4);
        alpha = 0.6;
    elseif nargin == 6
        alpha = 0.6;
    end
    
    
    %defines the randomness of searching strategy
    expT = 1.0e-5;
    
    %% learning algorithm
    maxIt = 1e1;
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
    currentQuality = eval(caseRs, currentR.*currentR, caseAlleleFreq, currentAlleleFreq, alpha);
    
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
        newQuality = eval(caseRs, newRs, caseAlleleFreq, newP, alpha);
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

%evaluate normalized r square difference and p difference 
function [newQuality] = eval(targetRs, currentRs, targetFreq, currentFreq, alpha)
    %we think R square has more weight than single allele frequence
    if nargin == 4
        alpha = 0.6;
    end

    pDiff = abs(targetFreq - currentFreq);
    normalDiff = sum(pDiff)/length(targetFreq);   %normalize
    
    rDiff = calcRDiff(targetRs, currentRs);
    [m, n] = size(targetRs);
    if m ~= n
        e = MException('blockReconstruct:x', 'in consistent dimenstion');
        throw(e);
    end
    nElements = m*(m-1)/2;
    normalRDiff = rDiff/nElements;
    
    newQuality = alpha*normalRDiff + (1-alpha)*normalDiff;
    newQuality = 100*newQuality;
end

function [newQuality] = calcRDiff(targetRs, newRs, smallRFilter)
    
    if nargin == 2
        smallRFilter = 0.0;
    end
    %make sure remove the diagnal elements
    targetRs(logical(eye(size(targetRs)))) = 0;
    newRs(logical(eye(size(newRs)))) = 0;
    
    %filter the small r values
    targetRs(targetRs < smallRFilter) = 0;
    newRs(targetRs < smallRFilter) = 0;
    
    newQuality = sum(sum(abs(targetRs - newRs)))/2;
end