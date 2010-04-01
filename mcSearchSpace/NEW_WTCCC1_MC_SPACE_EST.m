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
pseudoCount = 1;
%end configuration

cd(conf.dataPath);

rawFastaData = fastaread(conf.fastaFile);
[hap01Seq, alleleMapping] = encodeRawFastaSeq(rawFastaData);
[nInvariantSnps, indexs] = findInvariantSnps(hap01Seq);
[r c00 c01 c10 c11] = calcRfrom01seq(hap01Seq);

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
allbinslog = [bin0;bin1];
save('McEst.mat');

%% do analysis
allbins = [pow2(allbinslog(:,1)), allbinslog(:,2)];
[s idx] = sort(allbins(:,1), 'descend');
%the probability distribution
fp = allbins(idx,:);
cfp = cumsum(fp(:,1).*fp(:,2));
plot(cfp);
plot(fp(:,1).*fp(:,2), 'r.');

%find effective size cut
cutIdx = find(cfp>0.99999, 1, 'first');
effectiveSize = sum(fp(1:cutIdx,2));
fprintf('effective size = pow2(%d)\n', log2(effectiveSize));




