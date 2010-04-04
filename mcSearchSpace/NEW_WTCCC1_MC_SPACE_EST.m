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
conf.blocks = [1, 14;15, 25;26, 51;52, 78;79,119;120,136;137,154;155,177;178, 204];
conf.verbose = 1;
sampleSize = 100;%sequence to sample
acurateSteps = 12;
similarity = 0.9;
%---begin config---

%more bins, more acurate estimation
nBins = 2^acurateSteps;%distribution of probability, number of bins 
pseudoCount = 1;
%end configuration

cd(conf.dataPath);

rawFastaData = fastaread(conf.fastaFile);
[hap01Seq, alleleMapping] = encodeRawFastaSeq(rawFastaData);
[nInvariantSnps, indexs] = findInvariantSnps(hap01Seq);
r = calcRHapSeq(hap01Seq);
plotLD(r);

%profile on;
%build MC model from real fasta data
iMCmodel = iMCmodelBuild_N_state(hap01Seq, pseudoCount);
save('iMCmodel.mat', 'iMCmodel');
% try
%     %load('iMCmodel.mat');
% catch exception
%     iMCmodel = iMCmodelBuild_N_state(hap01Seq, pseudoCount);
%     save('iMCmodel.mat', 'iMCmodel');
% end

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
allbinslog = [bin0;bin1];
save('McEst.mat');


%%
% using block estimation
blocks = conf.blocks;
nBlock = size(conf.blocks, 1);
seqType = 1;
for i = 1:nBlock
    a = blocks(i, 1);
    b = blocks(i, 2);
    blockSeq = hap01Seq(:, a:b);
    u = unique(blockSeq, 'rows');
    nType = size(u, 1);
    seqType = seqType * nType;
end
blockEstSizelog = log2(seqType);
fprintf(1, 'block effectiveSize %f\n', blockEstSizelog);

% do analysis
allbins = [pow2(allbinslog(:,1)), allbinslog(:,2)];
[s idx] = sort(allbins(:,1), 'descend');
%the probability distribution
fp = allbins(idx,:);
cfp = cumsum(fp(:,1).*fp(:,2));

% find effective size and real size
cutIdx = find(cfp>0.99999, 1, 'first');
effectiveSize = sum(fp(1:cutIdx,2));
fprintf('effective size = pow2(%d)\n', log2(effectiveSize));
realSeqSize = 2^len;

% plot cfp
h = figure;
plot(cfp);
line([cutIdx, cutIdx], [0, 1.5], 'LineStyle', '--', 'Color', 'red');
title('cumulative distribution of p');
xlabel('bins');
ylabel('cumulative distribution of p');
legend('cdf', '0.9999 cut');
saveas(h, 'cfp.epsc');

%plot(fp(:,1).*fp(:,2), 'r.');

%find similarity space at alpha level
alphaSnps = floor((1-similarity)*len);
simSpace = 0;
for i = 1:alphaSnps
    simSpace = simSpace + nchoosek(len, i);
end
fprintf(1, 'simSpace = pow2(%f)', log2(simSpace));

%plot confliction index
K = 25;
avgConflictIndex = zeros(K, 3);
for k = 1:K
    sampleSize = k*len;
    mcSolutionlog = sampleSize*log2(effectiveSize) - log2(perms(sampleSize));
    realSolutionlog = sampleSize*len - log2(perms(sampleSize));
    blockSolutionlog = sampleSize*blockEstSizelog - log2(perms(sampleSize));
    constrainSpacelog = nchoosek(len, 2)*log2(sampleSize);
    avgConflictIndex(k, 1) = mcSolutionlog - constrainSpacelog;
    avgConflictIndex(k, 2) = realSolutionlog - constrainSpacelog;
    avgConflictIndex(k, 3) = blockSolutionlog - constrainSpacelog;
end

h = figure;
hold on;
plot(avgConflictIndex(:,1), 'rx-');
plot(avgConflictIndex(:,2), 'bo-');
plot(avgConflictIndex(:,3), 'g.-');
line([1, k], [0, 0]);
legend('MC Model', 'REAL', 'Block Est', 2);
title('avg unique solutions');
xlabel('sample size x*nSnps');
ylabel('confliction index');
hold off;
saveas(h, 'conflictIndex.epsc');