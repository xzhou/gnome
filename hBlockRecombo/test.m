%simpleLearning
cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77'

refFileName = 'SIM_100x77_ctl.fasta';
sampleFileName = 'SIM_100x77_smp.fasta';

debug = 1;

targetSeq = readSeq4(sampleFileName);
refSeq = readSeq4(refFileName);

alleleMapping = getMajorAllele(refSeq);

targetR = calcR(targetSeq, alleleMapping);
refR = calcR(refSeq, alleleMapping);

%use the major allele mapping

% b1 = [34 41];   %block one
% b2 = [42 45];   %block two haplotypeFreq = 

%     [hbtarget, ftarget] = getHaplotype(targetSeq, b1);
%     [hbref, fref] = getHaplotype(refSeq, b1);
%     
%     x = [hbtarget ftarget];
%     y = [hbref fref];
% 
%     x = sortrows(x);
%     y = sortrows(y);
    %b1 = [34 41];   %block one
    %b2 = [42 45];   %block two haplotypeFreq

    %b1 = [34 41];
    %b2 = [49 57];
    
    %b1 = [12 17];
    %b2 = [18 22];
    
    %b1 = [42 45];
    %b2 = [46 49];
result = zeros(100, 2);
allSvsQ = zeros(100, 2);

blocks = [12 18 34 42 46 49;
          17 22 41 45 49 57];

% for i = 1:100
%     [finalSeq, finalR, finalSignRate, SvsQ] = hbRecomb(targetR, refSeq, alleleMapping);
%     result(i, :) = [i finalSignRate];
%     allSvsQ(i,:) = SvsQ;
% end

finalSignRate = 0;
while finalSignRate < 0.8
    [finalSeq, finalR, finalSignRate, SvsQ] = hbRecomb(targetR, refSeq, alleleMapping, blocks(:, [1 2]));
end
maskMatrix1 = getBlockMask(targetR, blocks(:, [1 2]));
targetR(maskMatrix1) = targetR(maskMatrix1).*reshape(sign(finalR), [], 1);

finalSignRate = 0;
while finalSignRate < 0.8
    [finalSeq, finalR, finalSignRate, SvsQ] = hbRecomb(targetR, refSeq, alleleMapping, blocks(:, [3 4]));
end
maskMatrix2 = getBlockMask(targetR, blocks(:, [3 4]));
targetR(maskMatrix2) = targetR(maskMatrix2).*reshape(sign(finalR), [], 1);

finalSignRate = 0;
while finalSignRate < 0.8
    [finalSeq, finalR, finalSignRate, SvsQ] = hbRecomb(targetR, refSeq, alleleMapping, blocks(:, [4 5]));
end
maskMatrix3 = getBlockMask(targetR, blocks(:, [4 5]));
targetR(maskMatrix3) = targetR(maskMatrix3).*reshape(sign(finalR), [], 1);

maskMatrix = maskMatrix1 | maskMatrix2;
maskMatrix = maskMatrix | maskMatrix3;

% [finalSeq, finalR, finalSignRate, SvsQ] = hbRecomb(targetR, refSeq, alleleMapping, blocks(:, [3 6]));
% maskMatrix = getBlockMask(targetR, blocks(:, [3 6]));
% targetR(maskMatrix) = targetR(maskMatrix).*sign(finalR);

signMatrix = sign(targetR);

MyDrTest(1, 1, maskMatrix, signMatrix);

%hist(result(:,2), 0:0.1:1);
%save('resultrs.map', 'result');
%scatter(allSvsQ(:,1), allSvsQ(:,2));