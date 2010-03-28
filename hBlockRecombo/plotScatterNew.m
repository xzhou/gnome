function [ h ] = plotScatterNew(caseSeq4, refSeq4, testSeq4, caseR, refR, figTitle, iBigRepeat)
%in this function, we will plot the scatter graph based on 

if nargin <= 6
	figTitle = 'unknown title';
end

nS = length(caseSeq4(:,1));
Len = length(caseSeq4(1,:));

alleleMapping = getMajorAllele(refSeq4);

int2S = (caseSeq4 == repmat(alleleMapping,nS,1)) + 0;
int2R = (refSeq4 == repmat(alleleMapping,nS,1)) + 0;
int2T = (testSeq4 == repmat(alleleMapping,nS,1)) + 0;

index1 = [1: nS];
index2 = [nS+1: nS*2];
index3 = [nS*2+1:nS*3];

StatS.Tr = zeros(nS, 1);
StatR.Tr = zeros(nS, 1);
StatT.Tr = zeros(nS, 1);


for i = 1:nS
    StatS.Tr(i) = getTr(int2S(i,:), caseR, refR);
    StatR.Tr(i) = getTr(int2R(i,:), caseR, refR);
    StatT.Tr(i) = getTr(int2T(i,:), caseR, refR);
end
StatS.Tr = StatS.Tr/sqrt(Len*(Len-1)/2);
StatR.Tr = StatR.Tr/sqrt(Len*(Len-1)/2);
StatT.Tr = StatT.Tr/sqrt(Len*(Len-1)/2);

sortT = sort(StatT.Tr);

h = figure;
hold on;

plot(index1, StatS.Tr, '.r');
plot(index2, StatR.Tr, '.g');
plot(index3, StatT.Tr, '.b');
legend({'case' 'ref' 'test'});

plot(ones(3*nS)*sortT(int16(nS*0.95)));
plot(ones(3*nS)*sortT(int16(nS*0.99)));
xlabel('individual index');
ylabel('T_r value');
title(figTitle);
filename = strcat('Trial', num2str(iBigRepeat),figTitle);
print(h,'-dpdf',filename);

end

