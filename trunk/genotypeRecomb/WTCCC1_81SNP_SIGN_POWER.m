function [signPower]  = WTCCC1_81SNP_SIGN_POWER()
change_env();
startParallel();
dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/TPED';
fastaFile = 'Affx_gt_58C_Chiamo_07.tped.fasta';
alpha = 0.05;

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

caseSize = 100;
refSize = 100;
trial = 1000;
[caseTr, refTr] = sign_power(hap01Seq, caseSize, refSize, trial);

save('sign_power.mat');
load('sign_power.mat');

[nSeg, nTrial, nCase] = size(caseTr);
[~, ~, nRef] = size(refTr);

%0.95 identification rate
result = zeros(nSeg, nTrial);
for i = 1:nSeg
    %p = i * 0.1;
    for j = 1:nTrial
        caseTr_j = caseTr(i, j, :);
        refTr_j = refTr(i, j, :);
        id = floor((1-alpha)*nRef);
        refTr_j_sort = sort(refTr_j);
        T = refTr_j_sort(id);
        result(i, j) = sum(caseTr_j > T);
    end
end

maxIdentificationRate = max(result, [], 2);
avgIdr = mean(result, 2);
save('maxIdentificationRate.mat');

plot(0.1:0.1:1, [maxIdentificationRate, avgIdr]);
title(['signRate vs identification rate, nCase = ', num2str(nCase), ' alpha = ', num2str(alpha), ' tiral = ', num2str(trial)]);
xlabel('signRate')
ylabel('identification rate');
legend('max', 'mean');

%plot
% for i = 1:trial
%     p = i*0.1; 
%     MA = reshape(caseTr(i, 1, :), numel(caseTr(i, 1, :)), 1);
%     M0 = reshape(refTr(i, 1, :), numel(refTr(i, 1, :)), 1);
%     plot2hist(MA, M0, ['sign rate = ', num2str(p)], 10);
% end

end