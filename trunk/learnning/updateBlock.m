function d = updateBlock(A,r,fun,first,last)

% A: individual sequences
% r: rsquares for A
% first: the index of the first dirty SNP (the SNP which has been changed)
% last: the index of the last dirty SNP

[row,col]=size(A);
for i=1:col
    for k=first:last
        if i >= first && i <= last && i > k, continue; end
        if i == k, continue; end
        
        [x11,x12,x21,x22]=computeCounts(A,i,k);
        r(i,k) = feval(fun,x11,x12,x21,x22);
        r(k,i) = r(i,k);
    end
end
d = r;