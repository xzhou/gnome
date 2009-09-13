function [nextBlock] = getNextBlock(blocks, minBlockSize)
    %get Next block will return the row number of the next block for
    %learnning, we first seach the neighbour blocks with block size
    %larger than bsize
    
    if nargin == 1
        minBlockSize = 2;
    end
    
    %linear search the next largest neighbour block
    currentLargestBlockSize = 0;
    
    [m n] = size(blocks);
    for i = 1:(m-1)
        if blocks(i,4) == 1 && blocks(i+1,4) == 0 && blocks(i+1,3) > currentLargestBlockSize
            nextBlock = i + 1;
            currentLargestBlockSize = blocks(i+1, 3);
        elseif blocks(i,4) == 0 && blocks(i+1, 4) == 1 && blocks(i, 3) > currentLargestBlockSize
            nextBlock = i;
            currentLargestBlockSize = blocks(i, 3);
        end
    end
    
    if currentLargestBlockSize >= minBlockSize
        return;
    end
    
    %search for the next largest non-neighbour block
    currentLargestBlockSize = 0;
    for i = 1:(m-1)
        if blocks(i,4) == 0 && blocks(i,3) > currentLargestBlockSize
            nextBlock = i;
            currentLargestBlockSize = blocks(i,3);
        end
    end
    
    if currentLargestBlockSize >= minBlockSize
        return;
    end
    
    %we can not find any neighbour and non-neighbour block, return an
    %arbitary snps
    for i = 1:(m-1)
        if blocks(i,4) == 0
            nextBlock = i;
            return;
        end
    end
end