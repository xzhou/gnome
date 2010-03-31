function [T] = getThreshold(val, alpha)
%estimate the alpha level, correct checked!
beta = 1 - alpha;
sortVal = sort(val);
m = length(val);
id = round(m*beta);
T = mean([sortVal(id-1), sortVal(id), sortVal(id+1)]);
end