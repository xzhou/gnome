function [idr, T] = signRateIdrRandom(caseSeq, testSeq, fdr, level, caseR, refR, blocks)
%we fix the block and assume the block sign are the same, then randomly
%assign 
[m n] = size(caseR);
blockMask  = getAllBlockMask(caseR,blocks);
T = zeros(1, level);
idr=zeros(1, level);
for i = 1:level
    idr_i =  zeros(1, 10);
    Ti = idr_i;
    for k = 1:10
        p = i/level;
        maskP = getMaskP(n, p);
        mask = ((maskP+1) + blockMask)>0 + 0;
        caseMaskR = caseR.*mask;
        idr_i(k) = getIdr(caseSeq, testSeq, fdr, caseMaskR, refR);
        Ti(k) = (sum(sum(mask == 1)) - sum(diag(mask)))/2/nchoosek(n,2);
    end
    idr(i) = mean(idr_i);
    T(i) = mean(Ti);
end
end
