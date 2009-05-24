function rerr = computePairwiseFrequencyError(r1,r2,b,e)

% r1: rsquares
% r2: rsquares
% b: the index of the first SNP
% e: the index of the last SNP

rerr = 0;
for i=b:e-1
    for k=i+1:e
        %weight = 1-r1(i,k); % weight for the error
        %rerr = rerr + weight*abs(r1(i,k)-r2(i,k));
        rerr = rerr + abs(r1(i,k)-r2(i,k));
        %rerr = rerr + abs((r1(i,k)-r2(i,k)/r2(i,k))^2);
    end
end
len = e-b+1;
len = len * (len-1) / 2;
rerr = rerr / len;