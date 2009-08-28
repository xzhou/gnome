function [r] = calcR(int4seq, alleleMapping)
    %majorAlleleMapping = getMajorAllele(int4seq);
    %global majorAlleleMapping;
    %alleleMapping
    int2seq = (int4seq == repmat(alleleMapping, length(int4seq), 1)) + 0;
    r = corrcoef(int2seq);
    r(isnan(r)) = 0;
end
