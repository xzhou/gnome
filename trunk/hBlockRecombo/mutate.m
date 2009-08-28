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
