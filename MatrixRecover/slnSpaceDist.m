function slnSpaceDist(fileName)
% this function extract the solution distribution
close all;
if nargin == 0
    fileName = 'slnct_40-40:7-7x3000 2010-11-11 23:15:04.log';
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
            
%             outlierLim = 40;
%             filteredM = filteredM(filteredM(:,3)<outlierLim, :);
            
            maxSln = max(filteredM(:,3));
            bins = 0:20:maxSln;
            h = figure;
            counts = histc(filteredM(:,3), bins);
            %bar(counts(2:end));
            ratio = spaceRatioExact(m, n);
            hist(filteredM(:,3));

            %standard deviation
            s = std(filteredM(:,3));
            avg = mean(filteredM(:,3));
            titleStr = ['solution distribution ', 'M=', num2str(m), ' L=', ...
                num2str(n), ' repeat=', num2str(size(filteredM, 1)), ...
                ' ratio = ', num2str(ratio), ' avg =', num2str(avg), ' std = ', num2str(s)];
            
            vline(avg, 'r', 'mean');
            vline(avg-s, 'g', '\sigma');
            vline(avg-2*s, 'g', '\sigma');
            vline(ratio, 'b', 'est');
            
            title(titleStr);
            xlabel('number of solutions');
            ylabel('count');
            %ylim([0, size(filteredM,1)]);
            
            saveas(h, [titleStr, '.epsc']);
        end
    end
end