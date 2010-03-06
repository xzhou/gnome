function [c00] = calcC00(int4q, alleleMapping)
%   c00 is the major/major
[nIndividual, ~] = size(int4q);
    int2seq = (int4q == repmat(alleleMapping, nIndividual, 1)) + 0;
    [m n] = size(int2seq);
    c00 = ones(n, n);
    for i = 1:n-1
        for j = i:n
            snp1 = int2seq(:,i);
            snp2 = int2seq(:,j);
            ac00 = sum((snp1 + snp2) == 2);
            c00(i,j) = ac00;
        end
    end
    
    c00 = copyUpperToLower(c00);
    
end