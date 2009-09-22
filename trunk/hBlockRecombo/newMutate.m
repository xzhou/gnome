function [currentSeq] = newMutate(currentSeq, newBlock)
    [m n] = size(currentSeq);
    a = newBlock(1,1);  %start snps
    b = newBlock(1,2);  %end snps
    
    idx1 = randi([1 m]);
    idx2 = randi([1 m]);
    
    %mutate
    while currentSeq(idx1, a:b) == currentSeq(idx2, a:b)
        idx1 = randi([1 m]);
        idx2 = randi([1 m]);
    end
    
    %swap the sequence
    temp = currentSeq(idx1, a:b);
    currentSeq(idx1, a:b) = currentSeq(idx2, a:b);
    currentSeq(idx2, a:b) = temp;
end