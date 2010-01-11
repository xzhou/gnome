function [result] = getCounts(genotypeSeq)
    [m n] = size(genotypeSeq);

    pA = zeros(n,1);      %single allele frequency
    cA = zeros(n,1);
    
    % 9 genotype counts, eg. 
    % counts(1,1) 
    counts = zeros(n, n, 3, 3);
    
    %% count pA
    for i = 1:n-1
        for j = (i+1):n
            seqA = genotypeSeq(:,i);
            seqB = genotypeSeq(:,j);
            % some redundancy here
            cA(i) = 2*sum(seqA == 0) + sum(seqA == 1);
            cA(j) = (2*sum(seqB == 0) + sum(seqB == 1));
            pA(i) = (2*sum(seqA == 0) + sum(seqA == 1))/2.0/m;
            pA(j) = (2*sum(seqB == 0) + sum(seqB == 1))/2.0/m;
            
            
            % for pairwise
            for k = 1:m
                A = seqA(k);
                B = seqB(k);
                counts(i,j,A+1,B+1) = counts(i,j,A+1,B+1) + 1;
                counts(j,i,B+1,A+1) = counts(j,i,B+1,A+1) + 1;
            end
        end
    end
    
    result.pA = pA;
    result.counts = counts;
    result.cA = cA;
end