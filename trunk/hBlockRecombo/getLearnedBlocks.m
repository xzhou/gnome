function [nLearnedBlocks totalBlocks] = getLearnedBlocks(blocks)
    nLearnedBlocks = sum(blocks(:,4));
    [totalBlocks ncol] = size(blocks);
end