classdef GnomeCalculator
% GnomeCalculator will store the most basic method for gnome process    
    properties
    end
    
    methods
        
    end
    
    methods(Static)
        %get single allele frequency
        function [p, c] = getSingleAlleleFreq(sequence4, alleleMapping)
            sequence01 = GnomeCalculator.encode(sequence4, alleleMapping);
            [m,n] = size(sequence4);
            c = sum(sequence01);
            p = c*1.0/m;
        end
        
        %encode gnome of 4 int to sequence of 2
        function [sequence] = encode(seq, alleleMapping)
            [nIndividual, nSnps] = size(seq);
            sequence = (seq == repmat(alleleMapping, nIndividual, 1)) + 0;
        end
    end
end

