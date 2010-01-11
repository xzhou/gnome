classdef genotypeHelpFuncs
    %a wrapper of functoins for genotype learning
    properties
    end
    
    methods(Static)
        function [genotype] = readGenotypeFromFasta(fileName)
            seq4 = readSeq4(fileName);
            majorAllele = getMajorAllele(seq4);
            
            [m n] = size(seq4);
            if(mod(m, 2) ~= 0)
                disp 'unpaired sequence'
                m = m - 1;
                %e = MException('readGenotypeFromFasta:m', 'unpaired sequence');
                %throw(e)
            end
            
            nIndividual = round(m/2);
            
            genotype = zeros(nIndividual, n);
            for i = 1:nIndividual           
               for j = 1:n
                   a = seq4(2*i-1, j);%the first allele
                   b = seq4(2*i, j);%the next allele
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
    end
end

