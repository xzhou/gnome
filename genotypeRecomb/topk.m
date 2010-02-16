function [XYVal] = topk(M, K)
%sort and select the first K min num in matrix M
%XYVal = [x, y, value]

if nargin ~= 2 || numel(M) < K
    e = MException('topK:K', 'K not defined or numel is less thanK');
    throw(e);
end

XYVal = zeros(K, 3);
[m n] = size(M);

%to 1-D index
[linearM, id] = sort(reshape(M, 1, numel(M)));

linearM = linearM(1:K);
id = id(1:K);

%to 2-D index
j = floor((id-1)/m) + 1;
i = mod(id-1, m) + 1;

XYVal = [i', j', linearM'];
end

