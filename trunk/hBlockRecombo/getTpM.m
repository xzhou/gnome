function [tp] = getTpM(genoSeq, caseP, refP)
%calculate Homer's attacker's val, genoSeq must be encoded as in combineSeq
[m ~] = size(genoSeq);
tp = zeros(m,1);
for i = 1:m
   tp(i) = getTp(genoSeq(i,:), caseP, refP); 
end
end