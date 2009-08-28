function [seq4] = readSeq4(fileName)
    seq = fastaread(fileName);
    m = length(seq);
    n = length(seq(1).Sequence);
    seq4 = zeros(m, n);
    
    for i = 1:m
        seq4(i,:) = nt2int(seq(i).Sequence) - 1;
    end

end