function [r pA counts] = estimateR(genotypeSeq)
    [m n] = size(genotypeSeq);
    r = zeros(n, n);    %pairwise allele requency
    %allocate the data structure
    result = getCounts(genotypeSeq);
    
    pA = result.pA;
    counts = result.counts;
    
    for i = 1:n-1
       for j = i+1:n
            p1 = pA(i);
            p2 = pA(j);
            n3x3 = reshape(counts(i,j,:,:), 3, 3);
            x = mleR(p1, p2, n3x3);
            r(i,j) = x.r;
        end
    end
    
    %vectorize loop for speed
    for i = 1:n-1
        r(i+1:n,i) = r(i,i+1:n);
    end
end
