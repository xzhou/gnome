%using r square 
function [value sampleR] = evaluateSeq(targetRs, sampleSeq, alleleMapping)
    sampleR = calcR(sampleSeq, alleleMapping);
    sampleRs = sampleR.*sampleR;
    diff = abs(sampleRs - targetRs);
    diff(logical(eye(size(diff)))) = 0;
    value = sum(sum(diff.*diff))/2.0;
end