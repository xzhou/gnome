function [value, sampleR] = evaluateSeqAbs(targetRs, sampleSeq, alleleMapping)
    sampleR = calcR(sampleSeq, alleleMapping);
    targetR = sqrt(targetRs);
    diff = abs(abs(sampleR) - abs(targetR));
    diff(logical(eye(size(diff)))) = 0;
    value = sum(sum(abs(diff)))/2.0;
end