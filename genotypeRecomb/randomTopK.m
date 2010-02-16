function [xyVal] = randomTopK(M, K)
%function random top k using a random select according to the value in M,
%which is a matrix of numeric elements
%the frequency is inverse to it's value, less value, more likely to be
%selected

%1-D
M1 = reshape(M, 1, numel(M));
[SM1 id] = sort(M1);

nanFilter = ~isnan(SM1);
SM1 = SM1(nanFilter);
id = id(nanFilter);

Pr = 1./(sum(1./SM1).*SM1);
bins = cumsum(Pr);
xyVal = zeros(K, 3);
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
for i = 1:K
    r = rand();
    idx = whichBin(bins, r);
    [x y] = ind2sub(size(M), id(idx));
    xyVal(i,:) = [x y M(x, y)];
end
end

function [idx] = whichBin(bins, x)
%returns the index of bins that contains x
idx = -1;
bins = [0, bins];
m = numel(bins);
for i = 1:m
    if x < bins(i) && x >= bins(i-1)
        idx = i-1;
        return;
    end
end
end