function [blocks] = hbPartition(seq, alleleMapping, r)
%FINDHBLOCK find the the haplotype block from the r value
%We use dynamic programming algoirithm proposed by Zhang 2002
% to find the optimal block partition
    if nargin == 0
        error('not enough input');
    end
    
    %configuration
    delta = 0.0;
    
    m = length(seq);
    n = length(seq(1).Sequence);
    
    if nargin == 1
        %calculate r vlaue
        int4seq = zeros(m,n);
        
        parfor i = 1:m
            int4seq(i,:) = nt2int(seq(i).Sequence) - 1;
        end
        
        alleleMapping = getMajorAllele(int4seq);
        
        int2seq = (int4seq == repmat(alleleMapping, m, 1)) + 0;
        
        r = corrcoef(int2seq);
        
    end
        
    %cache for S function
    global smap partition intraCache interCache;
    smap = zeros(n, 2); %map index, value, partition pos
    partition = zeros(n, 1); %partition pos
    intraCache = zeros(n);
    interCache = zeros(n, n, n);
    
    S(r, n, delta);
    blocks = smap;
    
end

function [value, pos] = S(r, i, delta)
% S function calculate the partition for snps from 1:i
% r is the LD value with sign
% smap is the cache for dynamic programming
    global smap;
    if i == 0
        value = 0;
        pos = 0;
    elseif i == 1
        value = 0;
        pos = 0;
    elseif smap(i, 1) ~= 0
        value = smap(i, 1);
        pos = smap(i, 2);
    else
        %recursively calculate
        allj = zeros(i-1, 1);
        for j = 1:i - 1
            [jv jp] = S(r, j-1, delta);
            
            %we assume j is the new block
            intraj = intraBlockSum(r, j, i, delta);
            interj = interBlockSum(r, 1, j, i, delta);
            allj(j) = jv - interj + intraj;
        end
        [jv jpos] = max(allj(:,1));
        value = jv;
        pos = jpos;
        smap(i,:) = [value pos];
    end         
end

%calculate the sum of a internal block
function [value] = intraBlockSum(r, i, j, delta)
    global intraCache;
    if intraCache ~= 0
        value = intraCache(i,j);
    else
        subMatrix = r(i:j, i:j);
        subMatrix(logical(eye(j-i+1))) = 0;
        rsquare = subMatrix.*subMatrix - delta;
        value = sum(sum(rsquare.*rsquare))/2;
        intraCache(i,j) = value;
    end
end

%calculate the sum of inter block
function [value] = interBlockSum(r, low, mid, up, delta)
    global interCache;
    if interCache(low, mid, up) ~= 0
        value = interCache(low, mid, up);
    elseif low == mid
        value = 0;
    else
        leftSum = intraBlockSum(r, low, mid-1, delta);
        rightSum = intraBlockSum(r, mid, up, delta);
        total = intraBlockSum(r, low, up, delta);
        value = total - leftSum - rightSum;
        interCache(low, mid, up) = value;
    end
end