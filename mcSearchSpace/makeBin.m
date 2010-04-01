function [newBins, lowerBound] = makeBin(bins, nBins)
%create new bins that will do estimation of old bins to save sapce
%remove 

if any(bins(:,1) > 0)
    %all log probability should be less than 0
    e = MException('makeBin:p', 'log(p) <= 0');
    throw(e);
end

%count -Inf
cInf = sum(bins(bins(:,1)==-Inf, 2));
bins = bins(bins(:,1)~=-Inf, :);

lowerBound = floor(min(bins(:,1)));

binSize = abs(lowerBound/nBins);

%newBins = [avg_p, count]
newBins = zeros(nBins, 2);
%create bins
for i = 1:nBins
    binBtn = -i*binSize;
    binTop = binBtn + binSize;
    ps = bins(bins(:,1)>binBtn & bins(:,1)<=binTop, :);%HOT line
    if size(ps, 1) == 0
        newBins(i, :) = [binTop - binSize/2, 0];
    else
        counts = sum(ps(:,2));
        if(counts == 0)
            newBins(i,:) = [binTop - binSize/2, 0];
        else
            %avg_p = sum(ps(:,1).*ps(:,2))/counts;
            avg_p = log2(sum(pow2(ps(:,1)).*ps(:,2))/counts);
            newBins(i,:) = [avg_p, counts];
        end
    end
end

%newBins = [newBins; [-Inf, cInf]];

end