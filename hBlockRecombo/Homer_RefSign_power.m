%we need to know what's the power of copy sign and homer's test
if (~isdeployed)
    disp 'not deployed';
    cd('~/research_linux/gnome/bioWorkspace/genomeprj/common');
    change_env();
end

startParallel();
dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/fastPhase';
fastaFile = 'Affx_gt_58C_Chiamo_07.tped.200snp.extract.inp.fasta';
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

%% begin configuration
FDR = 0.05;
sameSize = 200;
caseSize = sameSize;
refSize = sameSize;
testSize = sameSize;
trial = 100;
nSnps = n;
useEstR = 1;

%% cut end snps
if nSnps < n
    nCut = n - nSnps;
    ub = n - floor(nCut/2);
    lb = floor(nCut/2) + 1;
    hap01Seq = hap01Seq(:, lb:ub);
end

%% try different case and reference
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));   %reset the random number generator

%idr = [Homer, Yong, Yong_copySign]
idr = zeros(trial, 3); %identification rate
%profile on;
parfor i = 1:trial
    %fprintf(1, 'trial %d\n', i);
    %randomly sample case, ref, test sample
    [caseSeq, refSeq, testSeq] = randomSampleCaseRefTest(hap01Seq, caseSize, refSize, testSize);
    caseGenoSeq = combineHapSeq(caseSeq);
    refGenoSeq = combineHapSeq(refSeq);
    testGenoSeq = combineHapSeq(refSeq);
    
    caseR = corrcoef(caseSeq);
    caseR(isnan(caseR)) = 0;
    refR = corrcoef(refSeq);
    refR(isnan(refR)) = 0;
    
    if useEstR == 1
        %calculate estimate R
        %fprintf(1, '\testR');
        caseEstR = estimateR(caseGenoSeq);
        caseEstR(isnan(caseEstR)) = 0;
        %fprintf(1, '.');
        refEstR = estimateR(refGenoSeq);
        refEstR(isnan(refEstR)) = 0;
        %fprintf(1, '.');
        %testEstR = estimateR(testGenoSeq);
        %fprintf(1, '.\n');
    else
        caseEstR = corrcoef(caseSeq);
        caseEstR(isnan(caseEstR)) = 0;
        refEstR = corrcoef(refSeq);
        refEstR(isnan(refEstR)) = 0;
    end
    
    %calculate single allele frequence
    caseP = sum(caseSeq)/caseSize/2;
    refP = sum(refSeq)/refSize/2;
    testP = sum(testSeq)/testSize/2;
    
    %calculate Homer's attack's identification rate
    Tp_case = getTpM(caseGenoSeq, caseP, refP);
    Tp_test = getTpM(testGenoSeq, caseP, refP);
    zp = getThreshold(Tp_test, FDR);%get percentile
    TpIdr = sum(Tp_case > zp);%get identification rate
    
    %calculate Yong's attack power
    Tr_case = getTrM(caseSeq, caseEstR, refEstR);
    Tr_ref = getTrM(refSeq, caseEstR, refEstR);
    Tr_test = getTrM(testSeq, caseEstR, refEstR);
    zr = getThreshold(Tr_test, FDR);%get fdr percentile
    TrIdr = sum(Tr_case > zr);%get identification rate

    %calculate Yong's attack using ref sign or copy sign
    caseEstRCopySign = caseEstR.*sign(refEstR);
    Tr_case_copySign = getTrM(caseSeq, caseEstRCopySign, refEstR);
    Tr_ref_copySign = getTrM(refSeq, caseEstRCopySign, refEstR);
    Tr_test_copySign = getTrM(testSeq, caseEstRCopySign, refEstR);
    zr_cs = getThreshold(Tr_test_copySign, FDR);
    TrIdr_cs = sum(Tr_case_copySign > zr_cs);
    
    %save result
    idr(i,:) = [TpIdr, TrIdr, TrIdr_cs];
    fprintf(1, '# %d\t%d\t%d\n', TpIdr, TrIdr, TrIdr_cs);
end
%profile viewer;
fileName = ['compare_power', num2str(caseSize), '.mat'];
save(fileName);

%%
%load(fileName);

%% plot result
%h = figure;
hold on;
h = plot(1:trial, idr(:,1), 'rx');
plot(trial+1:2*trial, idr(:,2), 'bo');
plot(2*trial+1:3*trial, idr(:,3), 'g.');
line([1,3*trial], [mean(idr(:,1)), mean(idr(:,1))], 'Color', 'red', 'LineWidth', 2);
line([1,3*trial], [mean(idr(:,2)), mean(idr(:,2))], 'Color', 'blue', 'LineWidth', 2);
line([1,3*trial], [mean(idr(:,3)), mean(idr(:,3))], 'Color', 'green', 'LineWidth', 2);
configStr = ['case', num2str(caseSize),'ref', num2str(refSize), ...
    'test', num2str(testSize), 'fdr', num2str(FDR), 'trial', num2str(trial), 'nSnps', num2str(nSnps), 'EstR', num2str(useEstR)];
title(configStr);
legend('Homer', 'All sign', 'Copy Sign', 2);
hold off;
saveas(h, [configStr, '.pdf']);


