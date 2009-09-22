function [finalSeq finalR finalSignRate] = SARecombo(targetR, targetSeq, currentSeq, blocks, newBlock, alleleMapping)
% function SARecombo do the combination

    rejHisotryLen = 10000;     %the iterations used to calculate the efficiency
    varLen = 5000;              %calculating the variance 
    initT = 10;               %the initial temperature
    maxRejections = 9900;        %the rejection
    
    coolingFactor = 0.9;
    k = 1.3806503E10-23;        %Boltzmann constant
    
    targetRs = targetR.*targetR;
    [initQuality currentR] = blockEvaluateSeq(targetRs, currentSeq, alleleMapping, blocks, newBlock);
    initSignRate = blockSignRate(targetR, currentR, blocks, newBlock);
    
    currentQuality = initQuality;
    signRate = initSignRate;
    T = initT;
    
    %record the rejection history
    rejHistory = zeros(rejHisotryLen);
    varHistory = zeros(varLen);
    
    rejHistoryIndex = 0;
    varHistoryIndex = 0;
    previousVar = 1;
    
    nRejection = 0;
    while nRejection <= maxRejections
        %record the rejection history
        rejHistoryIndex = mod(rejHistoryIndex+1, rejHisotryLen) + 1;
        varHistoryIndex = mod(varHistoryIndex+1, varLen) + 1;
        
        newSampleSeq = newMutate(currentSeq, newBlock);
        [newQuality newR] = blockEvaluateSeq(targetRs, newSampleSeq, alleleMapping, blocks, newBlock);
        delta = newQuality - currentQuality;
        if delta <= 0
            currentSeq = newSampleSeq;
            currentQuality = newQuality;
            currentR = newR;
            signRate = blockSignRate(targetR, newR, blocks, newBlock);
            rejHistory(rejHistoryIndex) = 0;
            varHistory(varHistoryIndex) = delta;
        else
            p = exp(-delta/k/T);
            x = rand(1);
            if x < p
                currentSeq = newSampleSeq;
                currentQuality = newQuality;
                currentR = newR;
                signRate = blockSignRate(targetR, newR, blocks, newBlock);
                rejHistory(rejHistoryIndex) = 0;
                varHistory(varHistoryIndex) = delta;
            else
                %rejected
                rejHistory(rejHistoryIndex) = 1;
            end
        end
        nRejection = sum(rejHistory);
        
        %TODO calculate the variance and cooling down
        
    end
    
    finalSeq = currentSeq;
    finalQual = currentQuality;
    finalR = currentR;
    finalSignRate = signRate;
    fprintf(1, 'initQ = %f \tfinalQ = %f initSR = %f \tfinalSR = %f\n', initQuality, finalQual, initSignRate, finalSignRate);
end