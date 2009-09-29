function [currentSeq] = shuffleNewBlock(currentSeq, block)
    [m n] = size(currentSeq);
    a = block(1, 1);
    b = block(1, 2);
    
    blockSeq = currentSeq(:,a:b);
    
    ordering = randperm(m);
    blockSeq = blockSeq(ordering, :);
    
    currentSeq(:,a:b) = blockSeq;
end