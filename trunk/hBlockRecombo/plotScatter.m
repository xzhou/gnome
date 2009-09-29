function [ h ] = plotScatter(caseSeq4, refSeq4, caseR, refR, figTitle)
%in this function, we will plot the scatter graph based on 

if nargin <= 4
	figTitle = 'unknown title';
end

nS = length(caseSeq4);
Len = length(caseSeq4(1,:));

alleleMapping = getMajorAllele(refSeq4);

int2S = (caseSeq4 == repmat(alleleMapping,nS,1)) + 0;
int2R = (refSeq4 == repmat(alleleMapping,nS,1)) + 0;

index1 = [1: nS];
index2 = [nS+1: nS*2];

StatS.Tr = zeros(nS, 1);
StatR.Tr = zeros(nS, 1);

for i = 1:nS
    StatS.Tr(i) = getTr(int2S(i,:), caseR, refR);
    StatR.Tr(i) = getTr(int2R(i,:), caseR, refR);
end
StatS.Tr = StatS.Tr/sqrt(Len*(Len-1)/2);
StatR.Tr = StatR.Tr/sqrt(Len*(Len-1)/2);

h = figure;
hold on;
plot(index1, StatS.Tr, '.r');
plot(index2, StatR.Tr, '.g');

xlabel('individual index');
ylabel('T_r value');
title(figTitle);
legend({'case' 'ref'});

end

