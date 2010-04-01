%zhou@indiana.edu Mar 16, 2010
%using MC model to estimate the feasible solution of matrix

%% begin configuration
OS = change_env();
startParallel();
if OS == 1
    conf.dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/fastPhase';
else
    %Yangyi, please add data path here
    conf.dataPath = '?';
end

%---begin config---
conf.fastaFile = 'Affx_gt_58C_Chiamo_07.tped.200snp.extract.inp.fasta';
conf.logFileName = 'WTCCC1_200_SNP_ANALYSIS.log';
conf.verbose = 1;
sampleSize = 100;%sequence to sample
acurateSteps = 12;

similarity = 0.9;
%---begin config---

%more bins, more acurate estimation
nBins = 2^acurateSteps;%distribution of probability, number of bins 
pseudoCount = 0;
%end configuration

cd(conf.dataPath);

rawFastaData = fastaread(conf.fastaFile);
[hap01Seq, alleleMapping] = encodeRawFastaSeq(rawFastaData);
[nInvariantSnps, indexs] = findInvariantSnps(hap01Seq);
%profile on;
%build MC model from real fasta data
try
    load('iMCmodel.mat');
catch exception
    iMCmodel = iMCmodelBuild_N_state(hap01Seq, pseudoCount);
    save('iMCmodel.mat', 'iMCmodel');
end

len = size(iMCmodel.transition, 3) + 1;
mcInit = iMCmodel.initial;
mcTransition = iMCmodel.transition;
nState = length(mcInit);


%the count of each probability bins
%bins = [probability, lastbits]
bins = zeros(nBins, 2);

%log all the probabilities to prevent zero-overflow
mcInit = log2(mcInit);
mcTransition = log2(mcTransition);


%start accurate steps
bin0 = mcInit(1);%init state 0
bin1 = mcInit(2);%init state 1

for i = 1:acurateSteps
    iTrans = mcTransition(:, :, i);
    newBin0 = [(bin0+iTrans(1, 1)); (bin1+iTrans(2,1))];
    newBin1 = [(bin0+iTrans(1, 2)); (bin1+iTrans(2,2))];
    bin0 = newBin0;
    bin1 = newBin1;
end

%start the estimation step
padOne = ones(length(bin0), 1);
bin0 = [bin0, padOne];%extend the probability with counts
bin1 = [bin1, padOne];

%re-bin
bin0 = makeBin(bin0, nBins);
bin1 = makeBin(bin1, nBins);
for i = acurateSteps+1:len-1
    fprintf(1, 'itreation = %d\n', i);
    iTrans = mcTransition(:, :, i);
    newBin0 = [[(bin0(:,1) + iTrans(1, 1)), bin0(:,2)]; [(bin1(:,1) + iTrans(2,1)), bin1(:,2)]];
    newBin1 = [[(bin0(:,1) + iTrans(1, 2)), bin0(:,2)]; [(bin1(:,1) + iTrans(2,2)), bin1(:,2)]];
    
    %re-bin
    bin0 = makeBin(newBin0, nBins);
    bin1 = makeBin(newBin1, nBins);
end

%% save old data
allbins = [bin0;bin1];
save('allbins.mat', 'allbins');

%% sort by probability
[s idx] = sort(allbins(:,1), 'descend');

%the distribution function of p
fp = allbins(idx, :);
logfp = log2(fp);
fc = fp(:, 2);
fccumsum = cumsum(fc);
x = fp(:,1);
h = plot(sampleSize*log2(fccumsum), 'b');
%set(gca, 'XDir', 'reverse');
hold on;
%plot real size
%plot([1, nBins], sampleSize*repmat(len, nBins, 1), 'r-.');
line([1;2*nBins], sampleSize*[len;len], 'Color', 'r');
%plot estimation
%plot([1, nBins], sampleSize*repmat(log2(fccumsum(end)), nBins, 1), 'g--');
line([1;2*nBins], sampleSize*[log2(fccumsum(end)), log2(fccumsum(end))], 'Color', 'g', 'LineStyle', '--');

title(['MC model effective searching space of WTCCC 1 CH7, sample size = ', num2str(sampleSize)]);
xlabel('bins');
ylabel('cumsum (log2)');
legend('cumsum', 'real size', 'effective size');
hold off;
saveas(h, 'effectiveSize.pdf');


%pdf of probability
figure;
ppdf = pow2(allbins(:,1)).*allbins(:,2);
ppdf = sort(ppdf, 'descend');
plot(cumsum(ppdf));
figure;
plot(ppdf);
save;
%calculate effective size

effectiveSeqSpace = fccumsum(end);
%snps that allowed to be the same

alphaSnp = floor((1-similarity)*len);
simSpace = 0;
for i = 1:alphaSnp
    simSpace = simSpace + nchoosek(len, i);
end

fprintf(1, 'simSpace = pow2(%f)', log2(simSpace));

%unique solutions at alpha level of similarity, estimation
K = 10;

avgConflictIndex = zeros(K, 2);
for k = 1:K
    sampleSize = k*len;
    mcSolutionslog = sampleSize*log2(effectiveSeqSpace) - log2(perms(sampleSize));
    realSolutionslog = sampleSize*len - log2(perms(sampleSize));
    constrainSpacelog = nchoosek(len, 2)*log2(sampleSize);
    avgConflictIndex(k, 1) = mcSolutionslog - constrainSpacelog;
    avgConflictIndex(k, 2) = realSolutionslog - constrainSpacelog;
end

h = plot(avgConflictIndex(:,1), 'rx-');
hold on;
plot(avgConflictIndex(:,2), 'bo-');
line([1, k], [0, 0]);
legend('MC Model', 'REAL');
title('avg uniquesolutions');
xlabel('sample size *nSnps');
ylabel('confliction index');
hold off;

fprintf(1, 'unique solutions = pow2(%f), realUnique = pow2(%f), similarity = %f\n', uniqueSolutionslog,realUniqueSolutionlog, similarity);



