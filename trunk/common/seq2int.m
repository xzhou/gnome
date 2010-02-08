function [seq4] = seq2int(fastaFormatSeq)
    m = length(fastaFormatSeq);
    n = length(fastaFormatSeq(1).Sequence);
    %seq4 = zeros(m, n);
    seq4 = struct('ID', [], 'Sequence', []);
    for i = 1:m
        seq4(i).ID = fastaFormatSeq(i).Header;
        seq4(i).Sequence = nt2int(fastaFormatSeq(i).Sequence) -1;
    end
end