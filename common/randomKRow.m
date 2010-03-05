function [subMatrix] = randomKRow(M, K)
%this function will randomly select K rows in M, no replacement
[m n] = size(M);
if m < K
    e = MException('randomKRow:error', 'not enough row');
    throw(e);
end
subMatrix = zeros(K, n);
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
a = m;
for i = 1:K
    id = randi(a, 1, 1);
    subMatrix(i,:) = M(id,:);
    a = a-1;
    M(id,:) = [];
end
end