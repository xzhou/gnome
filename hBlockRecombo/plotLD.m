function [h] = plotLD(r)
%plot LD structure using heat map
h = image(r.*r*100);
set(gca, 'YDir', 'normal');
load('MyColorMap', 'heatColorMap');
set(gcf, 'ColorMap', heatColorMap);
axis equal;
end