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
snpLen = 100;
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
hap01Seq = hap01Seq(:, 1:snpLen);
[nInvariantSnps, indexs] = findInvariantSnps(hap01Seq);
% r = calcRHapSeq(hap01Seq);
% plotLD(r);

%profile on;
%build MC model from real fasta data
iMCmodel = iMCmodelBuild_N_state(hap01Seq, pseudoCount);
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
save(['McEst', num2str(snpLen),'.mat']);


%%
% using block estimation
% blocks = conf.blocks;
% nBlock = size(conf.blocks, 1);
% seqType = 1;
% for i = 1:nBlock
%     a = blocks(i, 1);
%     b = blocks(i, 2);
%     blockSeq = hap01Seq(:, a:b);
%     u = unique(blockSeq, 'rows');
%     nType = size(u, 1);
%     seqType = seqType * nType;
% end
% blockEstSizelog = log2(seqType);
% fprintf(1, 'block effectiveSize %f\n', blockEstSizelog);
fsize = 22;

% do analysis
allbins = [pow2(allbinslog(:,1)), allbinslog(:,2)];
[s idx] = sort(allbins(:,1), 'descend');
%the probability distribution
fp = allbins(idx,:);
cfp = cumsum(fp(:,1).*fp(:,2));
cfp = cfp(find(cfp, 1):end);

% find effective size and real size
cutIdx = find(cfp>0.999999, 1, 'first');
effectiveSize = sum(fp(1:cutIdx,2));
fprintf('effective size = pow2(%d)\n', log2(effectiveSize));
realSeqSize = 2^len;

% plot cfp
h = figure;
set(gca, 'FontSize', fsize);
plot(cfp);
%line([cutIdx, cutIdx], [0, 1.1], 'LineStyle', '--', 'Color', 'red');
vh = vline(cutIdx, 'r', ['cut index = ', num2str(cutIdx)]);
%set(vh, 'FontSize', fsize);
%title('cumulative distribution of p', 'FontSize', 16);
xlabel('bins','FontSize', fsize);
ylabel('cumulative distribution of p','FontSize', fsize);
ylim([0, 1.1]);
legend('cdf', '0.99999 cut', 'FontSize', fsize);
set(gca, 'FontSize', fsize);
saveas(h, 'cfp.epsc');
saveas(h, 'cfp.fig');

% plot counts
h = figure;
set(gca, 'FontSize', fsize);
plot(log2(cumsum(fp(:,2))));
%line([cutIdx, cutIdx], [0, 1.1], 'LineStyle', '--', 'Color', 'red');
vh = vline(cutIdx, 'r', ['cut index = ', num2str(cutIdx)]);
xlabel('bins','FontSize', fsize);
ylabel('sequence counts (log_2)','FontSize', fsize);
legend('cumulative counts', 'FontSize', fsize);
set(gca, 'FontSize', fsize);
saveas(h, 'cfp_count.epsc');
saveas(h, 'cfp_count.fig');


%plot confliction index
K = 50;
avgConflictIndex = zeros(K, 3); %[MC, REAL, BLOCK]
SS = 2^len;
for k = 1:K
    sampleSize = k*len;
    mcSolutionlog = log2_mchoosen(effectiveSize, sampleSize);
    mcSolutionlog10 = log2_mchoosen(effectiveSize, 0.1*sampleSize);
    mcSolutionlog50 = log2_mchoosen(effectiveSize, 0.5*sampleSize);
    realSolutionlog = log2_mchoosen(SS, sampleSize);
    realSolutionlog10 = log2_mchoosen(SS, 0.1*sampleSize);
    realSolutionlog50 = log2_mchoosen(SS, 0.5*sampleSize);
    
    %TODO: 
    %blockSolutionlog = sampleSize*blockEstSizelog - log2(perms(sampleSize));
    constrainSpacelog = (nchoosek(len, 2)+1)*log2(sampleSize+1);
    avgConflictIndex(k, 1) = mcSolutionlog - constrainSpacelog;
    avgConflictIndex(k, 2) = mcSolutionlog10 - constrainSpacelog;
    avgConflictIndex(k, 3) = mcSolutionlog50 - constrainSpacelog; 
    avgConflictIndex(k, 4) = realSolutionlog - constrainSpacelog;
    avgConflictIndex(k, 5) = realSolutionlog10 - constrainSpacelog;
    avgConflictIndex(k, 6) = realSolutionlog50 - constrainSpacelog;
    %avgConflictIndex(k, 3) = blockSolutionlog - constrainSpacelog;
end



h = figure;
set(gca,'FontSize',fsize);
hold on;
[a, b] = getN(len, 0.5);
plot(avgConflictIndex(:,1), 'rx-');
%plot(avgConflictIndex(:,2), 'r+-');
plot(avgConflictIndex(:,3), 'ro-');
plot(avgConflictIndex(:,4), 'bs-');
%plot(avgConflictIndex(:,5), 'bd-');
plot(avgConflictIndex(:,6), 'b.-');

vline(ceil(a/len), 'g', '1');
vline(ceil(b/len), 'r', '0.5');

line([1, k], [0, 0]);
%legend('MC', 'MC0.1', 'MC0.5', 'REAL', 'REAL0.1', 'REAL0.5', 2);
legend('MC', 'MC 0.5', 'REAL', 'REAL0.5', 2);
%title('avg unique solutions', 'FontSize', 16);
xlabel('sample size x*L', 'FontSize', fsize);
ylabel('solutions (log2)', 'FontSize', fsize);

hold off;
saveas(h, 'conflictIndex.epsc');
saveas(h, 'conflictIndex.fig');