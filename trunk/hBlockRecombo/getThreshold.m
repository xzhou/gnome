function [T] = getThreshold(val, alpha)
%estimate the alpha level
beta = 1 - alpha;
sortVal = sort(val);
m = length(val);
id = m*beta;
T = sortVal(id);
end