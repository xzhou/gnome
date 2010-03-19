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

conf.fastaFile = 'Affx_gt_58C_Chiamo_07.tped.200snp.extract.inp.fasta';
conf.logFileName = 'WTCCC1_200_SNP_ANALYSIS.log';
conf.verbose = 1;

sampleSize = 100;
acurateSteps = 10;
nBins = 2^acurateSteps;%distribution of probability, number of bins 
pseudoCount = 0;
%end configuration

cd(conf.dataPath);

rawFastaData = fastaread(conf.fastaFile);
[hap01Seq, alleleMapping] = encodeRawFastaSeq(rawFastaData);

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
mcInit = log(mcInit);
mcTransition = log(mcTransition);


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

allbins = [bin0;bin1];

counts = exp(allbins(:,1)).*allbins(:,2);
counts = sort(counts, 'descend');
plot(cumsum(counts));



