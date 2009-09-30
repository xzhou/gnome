function [uniqueBLock, freq] = getBlockFreq(currentSeq, newBlock)
    [m n] = size(currentSeq);
    
    a = newBlock(1, 1);
    b = newBlock(1, 2);
    
    block = currentSeq(:,a:b);
    
    [uniqueBLock I J]  = unique(block, 'rows');
	
	freq = accumarray(J, 1);
    
end