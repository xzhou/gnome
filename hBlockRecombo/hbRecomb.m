function []= hbRecomb(targetR, refSeq)
% function hbRecomb use haplotype recombination to 
    %partition = hbPartition(refSeq);
    %assume we know the block structure
    
    maxIT = 10000000;
    maxDiff = 0.1;
    expT = 0.0001; %stacastic algorithm parameter
    
    
    b1 = [34 41];   %block one
    b2 = [42 45];   %block two haplotypeFreq = 
    
    targetR = targetR([b1(1):b1(2) b2(1):b2(2)], [b1(1):b1(2) b2(1):b2(2)]);
    
    targetRs = targetR.*targetR;
    
    refR = calcR(refSeq);
    refR = refR([b1(1):b1(2) b2(1):b2(2)], [b1(1):b1(2) b2(1):b2(2)]);
    
    targetRinit = targetR;
    refRinit = refR;
    
    [h1 f1] = getHaplotype(refSeq, b1);
    [h2 f2] = getHaplotype(refSeq, b2);
    
    lenB1 = b1(2) - b1(1) + 1;
    
    global majorAlleleMapping;
    
    
    currentSeq = [refSeq(:,b1(1):b1(2)) refSeq(:,b2(1):b2(2))];
    
    majorAlleleMapping = getMajorAllele(currentSeq);
    
    [currentQuality currentR] = evaluateSeq(targetRs, currentSeq);
    
    %sa algorithm
    disp('start stacastic learning');
    
    t = 0;
    while t < maxIT && currentQuality > maxDiff
        t = t + 1;
        newSampleSeq = mutate(currentSeq, lenB1);
        [newQuality newR] = evaluateSeq(targetRs, newSampleSeq);
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
                xdk = 0.0;
            end
        else
            %disp(t)
        end
    end
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


function [value] = evaluationR(targetRs, sampleRs)
    diff = abs(sampleRs - targetRs);
    diff(logical(eye(size(diff)))) = 0;
    value = sum(sum(diff.*diff))/2.0;
end

function [value sampleR] = evaluateSeq(targetRs, sampleSeq)
    sampleR = calcR(sampleSeq);
    sampleRs = sampleR.*sampleR;
    
    diff = abs(sampleRs - targetRs);
    diff(logical(eye(size(diff)))) = 0;
    value = sum(sum(diff.*diff))/2.0;
end


