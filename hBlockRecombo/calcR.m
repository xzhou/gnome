function [r] = calcR(int4seq, alleleMapping)
    %majorAlleleMapping = getMajorAllele(int4seq);
    %global majorAlleleMapping;
    %alleleMapping
    if nargin == 1
        alleleMapping = getMajorAllele(int4seq);
    end
    int2seq = (int4seq == repmat(alleleMapping, length(int4seq(:,1)), 1)) + 0;
    r = corrcoef(int2seq);
    r(isnan(r)) = 0;
end
