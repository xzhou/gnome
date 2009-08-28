function []= hbRecomb(targetR, refSeq, alleleMapping)
% function hbRecomb use haplotype recombination to infer the sign
    %partition = hbPartition(refSeq);
    %assume we know the block structure
    
    if nargin == alleleMapping
       alleleMapping = getMajorAllele(refSeq); 
    end
    
    maxIT = 10000;
    maxDiff = 0.1;
    expT = 0.005; %stacastic algorithm parameter
    
    %simulated 
    b1 = [34 41];   %block one
    b2 = [42 45];   %block two haplotypeFreq = 
    
    blockAlleleMapping = alleleMapping([b1(1):b1(2) b2(1):b2(2)]);
    
    refBlockSeq = refSeq(:,[b1(1):b1(2) b2(1):b2(2)]);
    refR = calcR(refBlockSeq, blockAlleleMapping);
    
    targetR = targetR([b1(1):b1(2) b2(1):b2(2)], [b1(1):b1(2) b2(1):b2(2)]);
    targetRs = targetR.*targetR;
    targetRinit = refR;
    refRinit = refR;
    
    
    [h1 f1] = getHaplotype(refSeq, b1);
    [h2 f2] = getHaplotype(refSeq, b2);
    
    lenB1 = b1(2) - b1(1) + 1;
    
    currentSeq = [refSeq(:,b1(1):b1(2)) refSeq(:,b2(1):b2(2))];
    
    [currentQuality currentR] = evaluateSeq(targetRs, currentSeq, blockAlleleMapping);
    initSignRate = SignRate(targetR, currentR);
    
    plotWithSignSimple(targetR, currentR);
    %sa algorithm
    disp('start stacastic learning');
    fprintf(1, '0\tinitQ = %f\t intSignRate = %f\n', currentQuality, initSignRate);
    
    t = 0;
    while t < maxIT && currentQuality > maxDiff
        t = t + 1;
        newSampleSeq = mutate(currentSeq, lenB1);
        [newQuality newR] = evaluateSeq(targetRs, newSampleSeq, blockAlleleMapping);
        Qdiff = newQuality - currentQuality;
        p = 1/(1+exp(Qdiff/expT));
        x = rand(1);
        if x < p
            currentSeq = newSampleSeq;
            currentQuality = newQuality;
            currentR = newR;
            signRate = SignRate(targetR, newR);
            fprintf(1, '%d \tQuality = %f\t SignRate = %f\n', t, currentQuality, signRate);
            if t >= 14500
                %a debug statementstatement
                xdebug = 0;
                pause;
            elseif signRate == 1.0
                %pause;
            end
        else
            %disp(t)
        end
    end
    
    plotWithSignSimple(targetR, currentR);

end

function [newSeq] = mutate(sampleSeq, lenB1)
    seq1 = sampleSeq(:, 1:lenB1);
    seq2 = sampleSeq(:, (lenB1+1):end);
    
    nSample = length(sampleSeq);
    
    idx1 = randi([1 nSample]);
    idx2 = randi([1 nSample]);
  
    while seq2(idx1,:) == seq2(idx2,:)
        idx1 = randi([1 nSample]);
        idx2 = randi([1 nSample]);
    end
    
    %swap
    temp = seq2(idx1,:);
    seq2(idx1,:) = seq2(idx2,:);
    seq2(idx2,:) = temp; 
    
    newSeq = [seq1 seq2];
end


function [value] = evaluationR(targetRs, sampleRs, alleleMapping)
    diff = abs(sampleRs - targetRs);
    diff(logical(eye(size(diff)))) = 0;
    value = sum(sum(diff.*diff))/2.0;
end

function [value sampleR] = evaluateSeq(targetRs, sampleSeq, alleleMapping)
    sampleR = calcR(sampleSeq, alleleMapping);
    sampleRs = sampleR.*sampleR;
    diff = abs(sampleRs - targetRs);
    diff(logical(eye(size(diff)))) = 0;
    value = sum(sum(diff.*diff))/2.0;
end


