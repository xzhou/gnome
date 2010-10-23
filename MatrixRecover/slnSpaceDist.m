function slnSpaceDist(fileName)
% this function extract the solution distribution
close all;
if nargin == 0
    fileName = 'slnct_10-11:3-40x300 2010-10-14 10:58:56.log';
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
        for j = 1:nn
            n = nSample(j);
            filteredM = M(M(:,2) == n, :);  %filtered m and n
            h = figure;
            hist(filteredM(:,3));
            
            ratio = spaceRatio(m, n);
            
            titleStr = ['solution distribution ', 'M=', num2str(m), ' L=', num2str(n), ' repeat=', num2str(size(filteredM, 1)), ' ratio = ', num2str(ratio)];
            title(titleStr);
            xlabel('number of solutions');
            ylabel('count');
            ylim([0, size(filteredM,1)]);
            saveas(h, [titleStr, '.pdf']);
        end
    end
end