function [nInvariant, idx] = findInvariantSnps(hap01seq)
[m n] = size(hap01seq);
index = (sum(hap01seq) == m);
x = 1:n;
idx = x(index);
nInvariant = length(idx);
end