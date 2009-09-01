function [haplotype, freq] = getHaplotype(seq, block)
% get haplotype block 
    l = block(1);
    u = block(2);
    blockSeq = seq(:, l:u);
    
    haplotype = unique(blockSeq, 'rows');
    
    [m n] = size(haplotype);
    freq = zeros(m, 1);

    for i = 1:m
        for j = 1:length(blockSeq(:,1))
            if blockSeq(j,:) == haplotype(i,:)
                freq(i,1) = freq(i,1) + 1;
            end
        end
    end
end