function slnSpaceDist(fileName)
% this function extract the solution distribution
if nargin == 0
    fileName = 'slnct_9-10:9-40x50.log';
end
    slnData = importdata(fileName);
    slnData = slnData.data;
    
%   [nInd, nSnp, probability]
    pMatrix = [];
    
    len = size(slnData, 1);
    
    mSample = unique(slnData(:,1));
    
    nm = size(mSample);
    
    for i = 1:nm
        m = mSample(i);
        M = slnData(slnData(:,1)==m, :);
        nSample = unique(M(:,2));
        nn = size(nSample);
        pMatrixOfm = [];
        for j = 1:nn
            n = nSample(j);
            filteredM = M(M(:,2) == n, :);  %filtered m and n
            hist(filteredM(:,3));
        end
    end
    hold off;
end