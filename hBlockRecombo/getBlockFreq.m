function [blockFreq] = getBlockFreq(currentSeq, newBlock)
%this function will count the number of different rows for a matrix and
%append the frequence at the last column of the different rows
    [m n] = size(currentSeq);
    
    a = newBlock(1, 1);
    b = newBlock(1, 2);
    
    block = currentSeq(:,a:b);
    
    [uniqueBlock I J]  = unique(block, 'rows');
	
	freq = accumarray(J, 1);
    
	blockFreq = [uniqueBlock freq];
end