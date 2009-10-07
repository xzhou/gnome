function [r] = calcBlockR(blocks, blockFreq, alleleMapping)
%function calcBlockR calculate the frequency of r from the block frequency
    seq = blockReconstruct(blocks, blockFreq, alleleMapping);
    r = calcR(seq, alleleMapping);
end

