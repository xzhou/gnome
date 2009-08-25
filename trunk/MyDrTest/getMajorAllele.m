function [majorAllele] = getMajorAllele(int4S)
    %returns the major allele frequency
    [m n] = size(int4S);
    counts = zeros(4,n);
    for i = 1:4
        counts(i,:) = sum(int4S == i-1);
    end
    
    maxCounts = max(counts);
    majorAllele = zeros(1,n);
    for i = 1:4
        index = counts(i,:) == maxCounts;
        majorAllele(index) = i - 1;
    end
end
