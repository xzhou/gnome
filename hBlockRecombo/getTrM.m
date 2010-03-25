function [tr] = getTrM(hapSeq, caseR, refR)
[m ~] = size(hapSeq);
tr = zeros(m, 1);
for i = 1:m
    tr(i) = getTr(hapSeq(i,:), caseR, refR);
end
end