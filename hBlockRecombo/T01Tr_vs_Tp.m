%this script will test the power of Tr and Homer's test over the simulated
%data

cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77';

allSeq = fastaread('80SNP_CEU_sim_4000seq_control1000.fasta');
caseSeq = readSeq4('SIM_100x77_ctl.fasta');
refSeq = readSeq4('SIM_100x77_smp.fasta');


T = 0;

% powerCurve(caseSeq, T, 200);

%use mask 
blocks = [12 18 34 42 46 49;
          17 22 41 45 49 57];
      
targetR = calcR(caseSeq);
refR = calcR(refSeq);

plotWithSignSimple(targetR, refR, 0.1);

signDiff = sign(targetR).*sign(refR);
signDiff(abs(targetR) <= 0.1) = 0;
totalSignDiff = sum(sum((signDiff == -1)+0))/2

[m n] = size(targetR);
[maskMatrix1 crossMatrix1] = getBlockMask(targetR, blocks(:, [1 2]));
[maskMatrix2 crossMatrix2] = getBlockMask(targetR, blocks(:, [3 4]));
[maskMatrix3 crossMatrix3] = getBlockMask(targetR, blocks(:, [4 5]));
[maskMatrix4 crossMatrix4] = getBlockMask(targetR, [50 55;54 57]);

crossMatrix = crossMatrix1 | crossMatrix2;
crossMatrix = crossMatrix | crossMatrix3;
crossMatrix = crossMatrix | crossMatrix4;

maskMatrix = maskMatrix1 | maskMatrix2;
maskMatrix = maskMatrix | maskMatrix3;
%maskMatrix = maskMatrix | maskMatrix4;

close all;

x = sum(sum(maskMatrix+0))/2
plotWithSignSimple(maskMatrix);

tempR = targetR;
tempR(tempR<0.1) = 0;
tempR(tempR>=0.1) = 1;
plotWithSignSimple(tempR);
nLE = (sum(sum(abs(targetR)>0.1))-n)/2

tempR = targetR;
idxMatrix1 = logical(triu(ones(size(targetR)), 6));
idxMatrix2 = logical(tril(ones(size(targetR)), -6));
tempR(idxMatrix1) = 0;
tempR(idxMatrix2) = 0;
tempR(tempR<0.1) = 0;
tempR(tempR>=0.1) = 1;
diagLE = ((sum(sum(abs(tempR)>=0.1)))-n)/2
plotWithSignSimple(tempR);

diagMask = ~(idxMatrix1 | idxMatrix2);

%all power
%powerCurve(allSeq, 0.0, 200, ones(size(targetR)));

%diag 0.1 non diag power
powerCurve(allSeq, 0.1, 200, diagMask);

