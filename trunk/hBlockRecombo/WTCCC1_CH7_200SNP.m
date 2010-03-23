%analysis WTCCC1 200 SNP and find the SNP structures

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

%end configuration

cd(conf.dataPath);
if conf.verbose
    conf.logfid = 1;
else
    conf.logfid = fopen(conf.logFileName, 'w');
end

%reading data
rawFastaData = fastaread(conf.fastaFile);
[hap01Seq, alleleMapping] = encodeRawFastaSeq(rawFastaData);
hapIntSeq = seq2int(rawFastaData);
hapIntSeqNoID = getSeqMatrix(hapIntSeq);
r = calcR(hapIntSeqNoID);
rs = r.*r;
%h = plotWithSignSimple(rs, rs);
h = plotWithSignSimple(abs(r), abs(r));
saveas(h, '200snp.LD.pdf');


