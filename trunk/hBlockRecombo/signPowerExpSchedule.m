%%experiments schedule

%% reading data
if (~isdeployed)
    disp 'not deployed';
    cd('~/research_linux/gnome/bioWorkspace/genomeprj/common');
    change_env();
end
startParallel();

%% data source
%dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/fastPhase';
%fastaFile = 'Affx_gt_58C_Chiamo_07.tped.200snp.extract.inp.fasta';

dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/fastPhase523';
fastaFile = 'Affx_gt_58C_Chiamo_07.tped.600SNP.extract.inp.fasta';

%% Yong's data source
% %real
% dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/YongExpData/real';
% fastaFile = 'chr10_FGFR2_200kb_phased_yri.fasta';
% 
% %sim
% dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/YongExpData/sim';
% fastaFile = '174SNP_CEU_sim_4000seq.fasta';

% dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/YongExpData/sim';
% fastaFile = '80SNP_YRI_sim_4000seq.fasta';
% 
% dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/YongExpData/sim';
% fastaFile = '80SNP_CEU_sim_4000seq.fasta';

dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/YongExpData/real';
fastaFile = 'hapmap_chr7_80SNP_CEU_haplotype.fasta';

%%init
logFile = 'signCompare.log';
logfid = fopen(logFile, 'w');
cd(dataPath);
fprintf(1, '%s\n', dataPath);

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

%% analysis LD structure
r = calcRHapSeq(hap01Seq);
save;
plotWithSignSimple(r.*r);

%% begin exp with schedule
fdrl = [0.05];
nSnps = [n];
sampleSize = [100, 200];%note this is individual size, sequence should 2*sampleSize
trials = 15;
levels = 10;
useEstR = 0;

for i = 1:length(fdrl)
    for j = 1:length(sampleSize);
        for k = 1:length(nSnps)
            powerAnalysis(hap01Seq, fdrl(i), sampleSize(j), sampleSize(j), sampleSize(j), trials, nSnps(k), levels, useEstR);
        end
    end
end