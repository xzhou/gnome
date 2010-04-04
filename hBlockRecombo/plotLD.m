function [h] = plotLD(r)
%plot LD structure using heat map
figure;
h = image(r.*r*100);
set(gca, 'YDir', 'normal');
%load('MyColorMap', 'heatColorMap');
heatColorMap = [0.972549021244049,0.972549021244049,0.972549021244049;0.975044548511505,0.975044548511505,0.884135484695435;0.977540135383606,0.977540135383606,0.795721948146820;0.980035662651062,0.980035662651062,0.707308351993561;0.982531189918518,0.982531189918518,0.618894815444946;0.985026717185974,0.985026717185974,0.530481278896332;0.987522304058075,0.987522304058075,0.442067742347717;0.990017831325531,0.990017831325531,0.353654175996780;0.992513358592987,0.992513358592987,0.265240639448166;0.995008885860443,0.995008885860443,0.176827087998390;0.997504472732544,0.997504472732544,0.0884135439991951;1,1,0;1,0.980769217014313,0;1,0.961538434028626,0;1,0.942307710647583,0;1,0.923076927661896,0;1,0.903846144676209,0;1,0.884615361690521,0;1,0.865384638309479,0;1,0.846153855323792,0;1,0.826923072338104,0;1,0.807692289352417,0;1,0.788461565971375,0;1,0.769230782985687,0;1,0.750000000000000,0;1,0.730769217014313,0;1,0.711538434028626,0;1,0.692307710647583,0;1,0.673076927661896,0;1,0.653846144676209,0;1,0.634615361690521,0;1,0.615384638309479,0;1,0.596153855323792,0;1,0.576923072338104,0;1,0.557692289352417,0;1,0.538461565971375,0;1,0.519230782985687,0;1,0.500000000000000,0;1,0.480769217014313,0;1,0.461538463830948,0;1,0.442307680845261,0;1,0.423076927661896,0;1,0.403846144676209,0;1,0.384615391492844,0;1,0.365384608507156,0;1,0.346153855323792,0;1,0.326923072338104,0;1,0.307692319154739,0;1,0.288461536169052,0;1,0.269230782985687,0;1,0.250000000000000,0;1,0.230769231915474,0;1,0.211538463830948,0;1,0.192307695746422,0;1,0.173076927661896,0;1,0.153846159577370,0;1,0.134615391492844,0;1,0.115384615957737,0;1,0.0961538478732109,0;1,0.0769230797886848,0;1,0.0576923079788685,0;1,0.0384615398943424,0;1,0.0192307699471712,0;1,0,0;];
set(gcf, 'ColorMap', heatColorMap);
axis equal;
end