%simpleLearning
cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77'

refFileName = 'SIM_100x77_ctl.fasta';
sampleFileName = 'SIM_100x77_smp.fasta';

debug = 1;

targetSeq = readSeq4(sampleFileName);
refSeq = readSeq4(refFileName);

alleleMapping = getMajorAllele(refSeq);

targetR = calcR(targetSeq, alleleMapping);
originalR = targetR;

refR = calcR(refSeq, alleleMapping);

case01 = targetR;
ref01 = refR;
ncaseLE01 = sum(sum(abs(case01)>0.1))/2;
nrefLE01 = sum(sum(abs(ref01)>0.1))/2;

signDiff = sign(case01.*ref01);
totaldiffRate = sum(sum(signDiff == -1))/2;

diffIndex = signDiff == -1;

%filter
case01Index = case01>=0.1;
ref01Index = ref01>=0.1;

x = case01Index&diffIndex + 0;

filterSignDiff = sum(sum(x));


result = zeros(100, 2);
allSvsQ = zeros(100, 2);

blocks = [12 18 34 42 46 49;
          17 22 41 45 49 57];

finalSignRate = 0;
while finalSignRate < 0.5
    [finalSeq, finalR, finalSignRate, SvsQ] = hbRecomb(targetR, refSeq, alleleMapping, blocks(:, [1 2]));
end
[maskMatrix1 crossMatrix1] = getBlockMask(targetR, blocks(:, [1 2]));
targetR(maskMatrix1) = abs(targetR(maskMatrix1)).*reshape(sign(finalR), [], 1);

finalSignRate = 0;
while finalSignRate < 0.8
    [finalSeq, finalR, finalSignRate, SvsQ] = hbRecomb(targetR, refSeq, alleleMapping, blocks(:, [3 4]));
end
[maskMatrix2 crossMatrix2] = getBlockMask(targetR, blocks(:, [3 4]));
targetR(maskMatrix2) = abs(targetR(maskMatrix2)).*reshape(sign(finalR), [], 1);

finalSignRate = 0;
while finalSignRate < 0.5
    [finalSeq, finalR, finalSignRate, SvsQ] = hbRecomb(targetR, refSeq, alleleMapping, blocks(:, [4 5]));
end
[maskMatrix3 crossMatrix3] = getBlockMask(targetR, blocks(:, [4 5]));
targetR(maskMatrix3) = abs(targetR(maskMatrix3)).*reshape(sign(finalR), [], 1);

maskMatrix = maskMatrix1 | maskMatrix2;
maskMatrix = maskMatrix | maskMatrix3;

crossMatrix = crossMatrix1 | crossMatrix2;
crossMatrix = crossMatrix | crossMatrix3;

% [finalSeq, finalR, finalSignRate, SvsQ] = hbRecomb(targetR, refSeq, alleleMapping, blocks(:, [3 6]));
% maskMatrix = getBlockMask(targetR, blocks(:, [3 6]));
% targetR(maskMatrix) = targetR(maskMatrix).*sign(finalR);

signMatrix = sign(targetR);

signDiff = sign(originalR).*sign(refR);
signDiff = signDiff.*maskMatrix;
signDiff(signDiff == 1) = 0;
signDiff = abs(signDiff);

%MyDrTest(1, 1, crossMatrix, signMatrix);
% powerTest(targetSeq, refSeq);       %the base test
% saveas(gcf, 'all_r_test.pdf');


%power test over 0.1
% maskMatrix = abs(originalR)>0.1;
% powerTest(targetSeq, refSeq, alleleMapping, targetR, maskMatrix, signDiff);
% 
% powerTest(targetSeq, refSeq, alleleMapping, targetR, maskMatrix, signMatrix);
% saveas(gcf, 'innerBlockWithCrossBlock.pdf');
% 
powerTest(targetSeq, refSeq, alleleMapping, targetR, crossMatrix, signMatrix);
saveas(gcf, 'crossBlcok.pdf');

powerTest(targetSeq, refSeq, alleleMapping, targetR, signDiff, signMatrix);
saveas(gcf, 'onlySignDiff.pdf');

%hist(result(:,2), 0:0.1:1);
%save('resultrs.map', 'result');
%scatter(allSvsQ(:,1), allSvsQ(:,2));