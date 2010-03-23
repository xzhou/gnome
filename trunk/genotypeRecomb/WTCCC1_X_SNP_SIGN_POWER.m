function [signPower]  = WTCCC1_X_SNP_SIGN_POWER()
%function used to estimate the sign agreement rate and power.
change_env();
startParallel();
dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/fastPhase';
fastaFile = 'Affx_gt_58C_Chiamo_07.tped.200snp.extract.inp.fasta';

alpha = 0.05;
cd(dataPath);

%reading fasta data
hapSeq = fastaread(fastaFile);
hapIntSeq = seq2int(hapSeq);
hapSeqNoID = getSeqMatrix(hapIntSeq);

alleleMapping = getMajorAllele(hapSeqNoID);

[m n] = size(hapSeqNoID);
hap01Seq = zeros(m, n);
for i= 1:m
    %1 major,  0 is minor
    hap01Seq(i,:) = (hapSeqNoID(i,:) == alleleMapping) + 0;
end

caseSize = 100;
refSize = 100;
nTest = 100;
trial = 100;
nSnps = 170;
useEstR = 1;
left = n - nSnps;
ub = n - left/2;
lb = left/2;

cuttedHap01Seq = hap01Seq(:, lb:ub);

[caseTr, refTr, testTr, caseTp, refTp, testTp] = sign_power(cuttedHap01Seq, caseSize, refSize, nTest, trial, 1);

save('sign_power.mat');
load('sign_power.mat');

[nSeg, nTrial, nCase] = size(caseTr);
[~, ~, nRef] = size(refTr);

%alpha FPR, Homer's attack, average
id = floor((1-alpha)*nTest);
sortTestTp = sort(testTp);
T = sortTestTp(id);
p_result = sum(caseTp > T);%base line

%0.95 FPR, we use test to find the .95 FPR line
result = zeros(nSeg, nTrial);
for i = 1:nSeg
    %p = i * 0.1;
    for j = 1:nTrial
        caseTr_j = caseTr(i, j, :);
        testTr_j = testTr(i, j, :);
        id = floor((1-alpha)*2*nTest);
        testTr_j_sort = sort(testTr_j);
        T = testTr_j_sort(id);
        result(i, j) = sum(caseTr_j > T);
    end
end

maxIdentificationRate = max(result, [], 2);
avgIdr = mean(result, 2);
save('maxIdentificationRate.mat');

h = plot(0.1:0.1:1, [maxIdentificationRate, avgIdr]);
hold on;
line([0, 1], [p_result, p_result]);
title(['signRate vs identification rate, nCase = ', num2str(nCase), ' nTest = ', num2str(nTest), ' alpha = ', num2str(alpha), ' tiral = ', num2str(trial), 'useEstR = ', num2str(useEstR)]);
xlabel('signRate')
ylabel('identification rate');
legend('max', 'mean');
hold off;
%saveas(h, 'signRate.pdf');


%plot 100% sign 
oneCaseTr = caseTr(10, 1, :);
oneTestTr = testTr(10, 1, :);
y = [reshape(oneCaseTr, 1, length(oneCaseTr)), reshape(oneTestTr, 1, length(oneTestTr))];
%scatter(1:length(y), y);

%plot
% for i = 1:trial
%     p = i*0.1; 
%     MA = reshape(caseTr(i, 1, :), numel(caseTr(i, 1, :)), 1);
%     M0 = reshape(refTr(i, 1, :), numel(refTr(i, 1, :)), 1);
%     plot2hist(MA, M0, ['sign rate = ', num2str(p)], 10);
% end

end