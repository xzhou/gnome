function [r_S single_S]= alleleMatch12(r_S, single_S, r_R, single_R)
% optimal matching of allel labels of two samples based on single SNP
% frequencies and r values.
% R is the reference
% S is the sample to be matched
% single.name: names for the two alles for each SNP
% single.p: frequenccies for the first allele
% r_S single_S will be adjusted based on the reference such that the
% simiarity is maximized

% decide the optimal matching based on single SNP frequencies
single_S.oldname = single_S.name;
single_S.name = single_R.name;
consistency = sign((single_R.p-0.5).*(single_S.p-0.5)); 
indx = find(consistency<0);

% for single SNP, modify S based on the matching
single_S.p(indx) = 1- single_S.p(indx);
for i = indx
    tmp = single_S.oldname{2,i};
    single_S.oldname{2,i} = single_S.oldname{1,i};
    single_S.oldname{1,i} = tmp;
end

% modify the sign of r_S values (LD values) based on the matching
consistency2 = consistency' * consistency; % whether the sign of r need to be changed for a pair
r_S = r_S.*consistency2;