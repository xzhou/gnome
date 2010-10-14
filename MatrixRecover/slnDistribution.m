function slnDistribution(fileName)
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
            filteredM = M(M(:,2) == n, :);
            nMultipleSln = sum(filteredM(:,3) > 1);
            nTotal = size(filteredM, 1);
            pMatrixOfm = [pMatrixOfm; [m, n, nMultipleSln*1.0/nTotal]];
        end
        pMatrix = [pMatrix; pMatrixOfm];
        
        hold on;
        plot(pMatrixOfm(:,2), pMatrixOfm(:,3));
    end
    hold off;
end