function [sequence] = encodeSequence(seq, alleleMapping)
%encode sequence will encode the standard ATCG sequence from haplotype 
%to 0/1
    nIndividual = length(seq);

    nSnps = length(seq(1).Sequence);
    
    int4S = zeros(nIndividual, nSnps);
    for i = 1:nIndividual
        int4S(i,:) = nt2int(seq(i).Sequence) - 1;
    end
    sequence = (int4S == repmat(alleleMapping, nIndividual, 1)) + 0;
end