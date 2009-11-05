function [r pA] = estimateR(genotypeSeq)
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
            n3x3 = counts(i,j,:,:);
            x = mleR(p1, p2, n3x3);
            r(i,j) = x.r;
        end
    end
end

function [result] = getCounts(genotypeSeq)
    [m n] = size(genotypeSeq);

    pA = zeros(n);      %single allele frequency
    
    % 9 genotype counts, eg. 
    % counts(1,1) 
    counts = zeros(n, n, 3, 3);
    
    %% count pA
    for i = 1:n-1
        for j = (i+1):n
            seqA = genotypeSeq(:,i);
            seqB = genotypeSeq(:,j);
            
            % some redundancy here
            pA(i) = (2*sum(seqA == 0) + sum(seqA == 1))/n;
            pA(j) = (2*sum(seqB == 0) + sum(seqB == 1))/n;
            
            % for pairwise
            for k = 1:n
                A = seqA(k);
                B = seqB(k);
                counts(i,j,A+1,B+1) = counts(i,j,A+1,B+1) + 1;
                counts(j,i,A+1,B+1) = counts(j,i,A+1,B+1) + 1;
            end
        end
    end
    
    result.pA = pA;
    result.counts = counts;
    
end