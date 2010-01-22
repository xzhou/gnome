function [value] = evaluationR(targetRs, sampleRs, alleleMapping)
    diff = abs(sampleRs - targetRs);
    diff(logical(eye(size(diff)))) = 0;
    value = nansum(nansum(diff.*diff))/2.0;
end