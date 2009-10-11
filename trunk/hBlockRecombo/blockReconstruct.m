function [seq] = blockReconstruct(blocks, blockFreq)
    [nBlock, snps] = size(blocks);
    [nFreq, tmp] = size(blockFreq);
    if nBlock ~= nFreq
        e = MException('blockReconstruct:x', 'in consistent dimenstion');
        throw(e);
    end
    
    rowIdentifier = zeros(sum(blockFreq),1);
    
    idx = 0;
    for i = 1:nFreq
       freq = blockFreq(i, 1);
       for j = 1:freq
           idx = idx + 1;
           rowIdentifier(idx,1) = i;
       end
    end
    seq = blocks(rowIdentifier,:);
    [m n] = size(seq);
    if m ~= sum(blockFreq)
        e = MException('blockReconstruct:x', 'can not recover the sequence!');
        throw(e);
    end
end