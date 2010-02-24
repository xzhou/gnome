function [r] = calcR(int4seq, alleleMapping, zeroNan)
    %majorAlleleMapping = getMajorAllele(int4seq);
    %global majorAlleleMapping;
    %alleleMapping
    if nargin == 1
        alleleMapping = getMajorAllele(int4seq);
    end
    if nargin < 3
        zeroNan = 0;
    end
    
    int2seq = (int4seq == repmat(alleleMapping, length(int4seq(:,1)), 1)) + 0;
    r = corrcoef(int2seq);
    
    if zeroNan == 0
        r(isnan(r)) = 0;
    end
end
