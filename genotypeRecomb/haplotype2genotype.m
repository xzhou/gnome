function [genotype] = haplotype2genotype(haplotypeSeq, majorAllele)
    if nargin == 1
        majorAllele = getMajorAllele(haplotypeSeq);
    end
    
    [m n] = size(haplotypeSeq);
    nIndividual = round(m/2);
    
    genotype = zeros(nIndividual, n);
    for i = 1:nIndividual           
       for j = 1:n
           a = haplotypeSeq(2*i-1, j);%the first allele
           b = haplotypeSeq(2*i, j);%the next allele
           m = majorAllele(j);%the major allele

           %convert to genotype
           if (a == m && b == m)
               genotype(i,j) = 0;
           elseif(a ~= m && b ~= m)
               genotype(i,j) = 2;
           else
               genotype(i,j) = 1;
           end
       end
    end
end