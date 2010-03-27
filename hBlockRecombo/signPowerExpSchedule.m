%%experiments schedule

%% reading data
if (~isdeployed)
    disp 'not deployed';
    cd('~/research_linux/gnome/bioWorkspace/genomeprj/common');
    change_env();
end

startParallel();
dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/fastPhase';
fastaFile = 'Affx_gt_58C_Chiamo_07.tped.200snp.extract.inp.fasta';

dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/fastPhase523';
fastaFile = 'Affx_gt_58C_Chiamo_07.tped.600SNP.extract.inp.fasta';

logFile = 'signCompare.log';
logfid = fopen(logFile, 'w');
cd(dataPath);

%% reading fasta data
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

fdrl = [0.01];
nSnps = [170, 200, 250, 300];
sampleSize = [200];
trials = 15;
levels = 10;

for i = 1:length(fdrl)
    for j = 1:length(sampleSize);
        for k = 1:length(nSnps)
            powerAnalysis(hap01Seq, fdrl(i), sampleSize(j), sampleSize(j), sampleSize(j), trials, nSnps(k), levels, 1);
        end
    end
end