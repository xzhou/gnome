function slnSpaceDist(fileName)
% this function extract the solution distribution
close all;
if nargin == 0
    fileName = 'slnct_10-10:5-6x30000 2010-10-31 01:14:14.log';
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
            maxSln = max(filteredM(:,3));
            bins = 0:1:maxSln;
            h = figure;
            counts = histc(filteredM(:,3), bins);
            bar(counts(2:end));
            ratio = spaceRatioExact(m, n);
            
            titleStr = ['solution distribution ', 'M=', num2str(m), ' L=', num2str(n), ' repeat=', num2str(size(filteredM, 1)), ' ratio = ', num2str(ratio)];
            title(titleStr);
            xlabel('number of solutions');
            ylabel('count');
            ylim([0, size(filteredM,1)]);
            saveas(h, [titleStr, '.pdf']);
        end
    end
end