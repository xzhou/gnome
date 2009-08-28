function [r] = calcR(int4seq)
    %alleleMapping = getMajorAllele(int4seq);
    global majorAlleleMapping;
    int2seq = (int4seq == repmat(majorAlleleMapping, length(int4seq), 1)) + 0;
    r = corrcoef(int2seq);
    r(isnan(r)) = 0;
end
