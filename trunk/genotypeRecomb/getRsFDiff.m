function [rsfDiff] = getRsFDiff(rs, f, targetRs, targetF, alpha)
%calculate the r square difference and f difference, after normalized to 1,
%the weight of f will be alpha
if nargin == 4
    alpha = 0.2;
end

[m n] = size(rs);

rsDiff = nansum(nansum(abs(diagZero(rs) - diagZero(targetRs))))/2;
rsDiffNorm = rsDiff*2/m/(m-1);

fDiff = sum(abs(f - targetF));
fDiffNorm = fDiff/m;
rsfDiff = (1-alpha)*rsDiffNorm + alpha*fDiffNorm;
end