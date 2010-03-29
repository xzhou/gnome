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
%filter1 =
%[39:62,94:112,123:142,166:179,190:202,203:210,249:256,259:277,299:313,323:331,450:478,480:497];
filter2 = [27:36,39:43,63:79,94:107,114:122,132:142,143:150,166:179,190:196,203:210,224:230, ...
240:248,259 277,289:297,315:319,323:331,345:353,450:478,480:497];

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

% dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/YongExpData/real';
% fastaFile = 'hapmap_chr7_80SNP_CEU_haplotype.fasta';

%%init
logFile = 'signCompare.log';
logfid = fopen(logFile, 'w');
cd(dataPath);
fprintf(1, '%s\n', dataPath);

%% reading fasta data
hapSeq = fastaread(fastaFile);
hapIntSeq = seq2int(hapSeq);
hapSeqNoID = getSeqMatrix(hapIntSeq);



%% filter
hapSeqNoID = hapSeqNoID(:, filter2);

%%
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
%h = plotLD(r);
%saveas(gcf, ['LD', num2str(n), '.fig']);

%% begin exp with schedule
fdrl = [0.01, 0.05];
nSnps = [n];
sampleSize = [200];%note this is individual size, sequence should 2*sampleSize
trials = 15;
levels = 10;
useEstR = 1;

for i = 1:length(fdrl)
    for j = 1:length(sampleSize);
        for k = 1:length(nSnps)
            powerAnalysis(hap01Seq, fdrl(i), sampleSize(j), sampleSize(j), sampleSize(j), trials, nSnps(k), levels, useEstR);
        end
    end
end