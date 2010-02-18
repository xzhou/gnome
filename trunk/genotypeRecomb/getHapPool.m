function [hapPool] = getHapPool(hapSeqs, blocks)
%extract the haplotype block structure from hap sequence, the block
%information is in config
[nSeq, nSnps] = size(hapSeqs);

if nargin == 1
    blocks = [1, nSnps];%defualt as one block
end

[nBlock, tmp] = size(blocks);

hapPool = struct('uniqueBlock',[], 'freq', []);

for i = 1:nBlock
    hapPool(i) = gBlockFreq(hapSeqs, blocks(i,:));
end

end