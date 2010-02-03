function [seq4] = seq2int(fastaFormatSeq)
    m = length(fastaFormatSeq);
    n = length(fastaFormatSeq(1).Sequence);
    seq4 = zeros(m, n);
    for i = 1:m
        seq4(i, :) = nt2int(fastaFormatSeq(i).Sequence) -1;
    end
end