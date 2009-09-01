function [maskMatrix] = getBlockMask(all_r, blocks)
    maskMatrix = zeros(size(all_r));
    a = blocks(1, 1);
    b = blocks(2, 1);
    c = blocks(1, 2);
    d = blocks(2, 2);
    
    maskMatrix(a:b, c:d) = 1;
    maskMatrix(c:d, a:b) = 1;
    maskMatrix(a:b, a:b) = 1;
    maskMatrix(c:d, c:d) = 1;
    maskMatrix = logical(maskMatrix);
end