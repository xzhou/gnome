function [genoSeq] = combineHapSeq(hapSeq)
%1 is major, 0 is minor
[nHapSeq, nSnps] = size(hapSeq);
genoSeq = zeros(nHapSeq/2, nSnps);
for i = 1:nHapSeq/2
    for j = 1:nSnps
        a1 = hapSeq(2*i,j);
        a2 = hapSeq(2*i-1, j);
        if a1 == 1 && a2 == 1
            genoSeq(i,j) = 0;
        elseif a1 ~= 1 && a2 ~= 1
            genoSeq(i,j) = 2;
        else
            genoSeq(i,j) = 1;
        end
    end
end
end
