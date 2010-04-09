function [idr, T, maxIdr, maxIdr2] = signRateIdrRandom(caseSeq, test2Seq, testSeq, fdr, level, caseR, refR, blocks)
%we fix the block and assume the block sign are the same, then randomly
%assign 
[m n] = size(caseR);
blockMask  = getAllBlockMask(caseR,blocks);
T = zeros(1, level);
idr=zeros(1, level);
maxIdr = idr;
maxIdr2 = idr;
for i = 1:level
    % idr_i =  zeros(1, 10);
    % idr_i2 =  zeros(1, 10);
    % Ti = idr_i;
	idr_i =  zeros(1, 1);
    idr_i2 =  zeros(1, 1);
    Ti = idr_i;
    for k = 1:1
        p = i/level;
        maskP = getMaskP(n, p);
        %keep the sign in the block 
        maskP(logical(blockMask)) = 1;
        caseMaskR = caseR.*maskP;
        idr_i(k) = getIdr(caseSeq, testSeq, fdr, caseMaskR, refR);
        idr_i2(k) = getIdr(test2Seq, testSeq, fdr, caseMaskR, refR);
        Ti(k) = (sum(sum(maskP == 1)) - sum(diag(maskP)))/2/nchoosek(n,2);
    end
    maxIdr(i) = max(idr_i);
	maxIdr2(i) = max(idr_i2);
    idr(i) = mean(idr_i);
    T(i) = mean(Ti);
end
end
