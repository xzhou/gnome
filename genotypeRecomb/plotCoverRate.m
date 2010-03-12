function [h] = plotCoverRate(freqInfo, tag)
%plot the cover rate and block type
if nargin == 1
    tag = '';
end

[nBlock ~] = size(freqInfo);

for i = 1:nBlock
    freq = freqInfo{i, 1}.freq;
    sortedFreq = sort(freq, 'descend');
    total = sum(freq);
    cumsumFreq = cumsum(sortedFreq/total);
    normalFreq = sortedFreq/total;
    h = plot(cumsumFreq);
    hold on;
    bar(normalFreq);
    name = [tag, 'block', num2str(i), '.pdf'];
    xlabel('number of genotype');
    ylabel('cover rate');
    title(['block', num2str(i)]);
    ylim([0, 1]);
    saveas(h, name);
    hold off;
end
end