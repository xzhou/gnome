function [targetR refSeq blkMapping] = getTarget(allTargetR, allRefSeq, allAlleleMapping, blocks)
%get the corresponding target R and reference Sequence
    a = blocks(1,1);
    b = blocks(1,2);
    
    c = blocks(2,1);
    d = blocks(2,2);
    
    
    %length of block 1
    bLen1 = b - a + 1;
    bLen2 = d - c + 1;
    
    b2Start = bLen1 + 1;
    b2End = bLen1 + bLen2;
    
    targetR = zeros(bLen1 + bLen2);
    
    %fill the diagnal data
    targetR(1:bLen1, 1:bLen1) = allTargetR(a:b, a:b);
    targetR(b2Start:b2End, b2Start:b2End) = allTargetR(c:d, c:d);
    
    %fill the cross data
    targetR(1:bLen1, b2Start:b2End) = allTargetR(a:b, c:d);
    targetR(b2Start:b2End, 1:bLen1) = allTargetR(c:d, a:b);
    
    %get the sequence
    refSeq = [allRefSeq(:, a:b) allRefSeq(:, c:d)];
    
    blkMapping = [allAlleleMapping(a:b) allAlleleMapping(c:d)];
end