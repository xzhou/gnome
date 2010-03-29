function [mask, p] = getAllBlockMask(r, blocks)
%%get the mask for all blocks and calculate the percentage of mask
[m ~] = size(blocks);
mask = zeros(size(r));
for i = 1:m
    a = blocks(i, 1);
    b = blocks(i, 2);
    mask(a:b, a:b) = 1;
end
diagSum = sum(blocks(:,2)-blocks(:,1));
p = (sum(sum(mask)) - diagSum)/(numel(r) - m);
end