function [ blockFreq ] = gBlockFreq(gSeq, newBlock)
% get the frequency of gnoetype blocks
    %[m n] = size(gSeq);
  a = newBlock(1, 1);
  b = newBlock(1, 2);
    
  block = gSeq(:,a:b);
    
  [uniqueBlock I J]  = unique(block, 'rows');
	
	freq = accumarray(J, 1);
    
	blockFreq.uniqueBlock = uniqueBlock;
  blockFreq.freq = freq;
end

