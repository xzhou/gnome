function [r c00 c01 c10 c11] = calcPairwiseFreq(int4q, alleleMapping)
%   c00 is the major/major
    int2seq = (int4q == repmat(alleleMapping, length(int4q), 1)) + 0;
    [m n] = size(int2seq);
    
    for i = 1:n
        for j = 1:n
            snp1 = int2seq(:,i);
            snp2 = int2seq(:,j);
            ac00 = sum((snp1 + snp2) == 2);
            ac01 = sum((snp1 - snp2) == 1);
            ac10 = sum((snp1 - snp2) == -1);
            ac11 = sum((snp1 + snp2) == 0);
            c00(i,j) = ac00;
            c01(i,j) = ac01;
            c10(i,j) = ac10;
            c11(i,j) = ac11;
            
            D = ac00*ac11 - ac01*ac10;
            L = (ac00 + ac01)*(ac00+ac10)*(ac11+ac01)*(ac11+ac10);
            
            if L == 0
                rij = 0;
            else
                rij = D/sqrt(L);
            end
            r(i,j) = rij;
        end
    end

end