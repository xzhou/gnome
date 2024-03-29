function [signPower]  = WTCCC1_81SNP_SIGN_POWER()
change_env();
startParallel();
dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/TPED';
fastaFile = 'Affx_gt_58C_Chiamo_07.tped.fasta';

cd(dataPath);

hapSeq = fastaread(fastaFile);
hapIntSeq = seq2int(hapSeq);
hapSeqNoID = getSeqMatrix(hapIntSeq);

alleleMapping = getMajorAllele(hapSeqNoID);

[m n] = size(hapSeqNoID);
hap01Seq = zeros(m, n);
for i= 1:m
    hap01Seq(i,:) = (hapSeqNoID(i,:) == alleleMapping) + 0;
end

caseSize = 300;
refSize = 100;
nTest = 100;
trial = 300;
alpha = 0.05;

[caseTr, refTr, testTr, ~, ~, ~, ~, ~] = sign_power(hap01Seq, caseSize, refSize, nTest, trial, 0);

save('sign_power.mat');
load('sign_power.mat');

[nSeg, nTrial, nCase] = size(caseTr);
[~, ~, nRef] = size(refTr);

%0.95 FPR, we use test to find the .95 FPR line
result = zeros(nSeg, nTrial);
for i = 1:nSeg
    %p = i * 0.1;
    for j = 1:nTrial
        caseTr_j = caseTr(i, j, :);
        testTr_j = testTr(i, j, :);
        id = floor((1-alpha)*nTest);
        testTr_j_sort = sort(testTr_j);
        T = testTr_j_sort(id);
        result(i, j) = sum(caseTr_j > T);
    end
end

maxIdentificationRate = max(result, [], 2);
avgIdr = mean(result, 2);
save('maxIdentificationRate.mat');

h = plot(0.1:0.1:1, [maxIdentificationRate, avgIdr]);
title(['signRate vs identification rate, nCase = ', num2str(nCase), ' nTest = ', num2str(nTest), ' alpha = ', num2str(alpha), ' tiral = ', num2str(trial)]);
xlabel('signRate')
ylabel('identification rate');
legend('max', 'mean');
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