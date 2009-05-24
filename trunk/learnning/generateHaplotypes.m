function hap=generateHaplotypes(iSNP)

% iSNP: the # of SNPs
% hap: haplotypes with iSNP SNPs
% For example, if iSNP is 3, hap will include following haplotypes:
% 0 0 0
% 0 0 1
% 0 1 0
% 0 1 1
% 1 0 0
% 1 0 1
% 1 1 0
% 1 1 1

hap=[];
dims=[];
for i=0:(2^iSNP-1)
    tmp = i;
    for k=(iSNP-1):-1:0
        dims(k+1) = fix(tmp/2^k);
        tmp = mod(tmp,2^k);
    end
    hap = [hap;dims];
end