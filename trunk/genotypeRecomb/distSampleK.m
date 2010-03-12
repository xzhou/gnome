function [val] = distSampleK(val, freq, K)
%% sample K values from val according to freq
n = sum(freq);
pool = generatePool(val, freq);
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
ix = randi(n, K, 1);
val = pool(ix,:);
end

function [seq] = generatePool(val, freq)
n = length(freq);
[~, mv] = size(val);
seq = zeros(sum(freq), mv);
count = 0;
for i = 1:n
    f = freq(i);
    a = count + 1;
    b = count + f;
    seq(a:b, :) = repmat(val(i,:), f, 1);
    count = count + f;
end
end