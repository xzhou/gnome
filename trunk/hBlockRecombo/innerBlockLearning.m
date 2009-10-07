function [] = innerBlockLearning(caseBlock, caseFreq, refBlock, refFreq, refcaseFreq, alleleMapping)
    if nargin == 5
        alleleMapping = getMajorAllele(refSeq4);
    end
    
    expT = 1.0e-10;
    
    %% learning algorithm
    maxIt = 100000;
    itr = 0;
    
    targetSeq = blockReconstruct(caseBlock, caseFreq);
    caseR = calcR(targetSeq, alleleMapping);
    caseRs = caseR.*caseR;
    
    currentSeq = blockReconstruct(refBlock, refFreq);
    currentR = calcR(currentSeq, alleleMapping);
    currentRs = currentR.*currentR;
    currentFreq = refFreq;
    currentQuality = eval(caseRs, currentR.*currentR);
    
    fprintf(1, '%f\n', currentQuality);
    
    while itr < maxIt
        itr = itr + 1;
        [newSeq4 newFreq] = getNextSeq(refBlock, currentFreq);
        newR = calcR(newSeq4, alleleMapping);
        newRs = newR.*newR;
        newQuality = eval(caseRs, newRs);
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
            fprintf(1, 'newQuality = %f\n', newQuality);
        else
            %do nothing
        end
    end
    
    %end of learning
    finalFreq = currentFreq;
    finalQual = currentQuality;
    finalR = currentR;
    
    save('innerBlockLearning.mat');
end

function [newSeq, newBlockFreq] = getNextSeq(hyplotypes, currentBlockFreq)
    newBlockFreq = blockNaiveMutate(currentBlockFreq);
    newSeq = blockReconstruct(hyplotypes, newBlockFreq);
end

function [newQuality] = eval(targetRs, currentRs)
    newQuality = calcRDiff(targetRs, currentRs);
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