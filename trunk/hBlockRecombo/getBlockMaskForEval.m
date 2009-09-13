function [blockMask] = getBlockMaskForEval(r,learnedBlocks, newBlock)
%GETBLOCKMASKFOREVAL 
%   @r: the frame of r value
    [m n] = size(learnedBlocks);
    
    blockMask = zeros(size(r));
    %for each pair of learned blocks
    for i = 1:m
        if learnedBlocks(i,4) ~= 1
            continue;
        end
        for j = (i+1):m
            if learnedBlocks(j,4) ~= 1
                continue;
            end
            
            start_i = learnedBlocks(i,1);   %start snps
            end_i = learnedBlocks(i,2);     %end snps
            start_j = learnedBlocks(j,1);   %2nd block start
            end_j = learnedBlocks(j,2);     %2nd block end
            
            %inner block
            blockMask(start_i:end_i, start_i:end_i) = 1;
            blockMask(start_j:end_j, start_j:end_j) = 1;
            
            %cross block mask
            blockMask(start_i:end_i, start_j:end_j) = 1;
            blockMask(start_j:end_j, start_i:end_i) = 1;      
        end
    end
    
    %add the new block
    if nargin == 3
        newBlockStart = newBlock(1);
        newBlockEnd = newBlock(2);
        for i = 1:m
            if learnedBlocks(i, 4) ~= 1
                continue
            end
            a = learnedBlocks(i, 1);
            b = learnedBlocks(i, 2);
            
            blockMask(a:b, a:b) = 1;
            blockMask(newBlockStart:newBlockEnd, newBlockStart:newBlockEnd) = 1;
            blockMask(newBlockStart:newBlockEnd, a:b) = 1;
            blockMask(a:b, newBlockStart:newBlockEnd) = 1;
        end
    end
    
end