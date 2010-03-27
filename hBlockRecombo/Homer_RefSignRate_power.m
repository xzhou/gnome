%we need to know what's the power of copy sign and homer's test
function [] = Homer_RefSignRate_power(FDR, nSnps, sampleSize)
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

%%config Pool
fdrl = [0.01, 0.05];
nSnps = [100, 200, 250, 300];
sampleSize = [200];

for i = 1:length(fdrl)
    for j = 1:length(sampleSize);
        for k = 1:length(nSnps)
            
        end
    end
end
%% begin configuration
FDR = 0.01;
sameSize = 200;
caseSize = sameSize;
refSize = sameSize;
testSize = sameSize;
trial = 20;
nSnps = 200;
useEstR = 1;
levels = 10; %divide sign recover rate by 10 level
configStr = ['case', num2str(caseSize),' ref', num2str(refSize), ...
    ' test', num2str(testSize), ' fdr', num2str(FDR), ' trial', num2str(trial), ' nSnps', num2str(nSnps), ' EstR', num2str(useEstR)];
fprintf(1, '%s', configStr);


%% cut end snps
if nSnps < n
%     nCut = n - nSnps;
%     ub = n - floor(nCut/2);
%     lb = floor(nCut/2) + 1;
%     hap01Seq = hap01Seq(:, lb:ub);
hap01Seq = hap01Seq(:, 1:nSnps);
end

%% try different case and reference
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));   %reset the random number generator

%idr = [signLevel, trials]
idrSignPower = zeros(levels, trial);
T = idrSignPower;
idrEstSignPower = zeros(levels, trial);
idrP = zeros(1, trial);
idrCS = zeros(1, trial);
idrCSEst = idrCS;

%profile on;
parfor i = 1:trial
    fprintf(1, 'trial %d\n', i);
    %randomly sample case, ref, test sample
    [caseSeq, refSeq, testSeq] = randomSampleCaseRefTest(hap01Seq, caseSize, refSize, testSize);
    caseGenoSeq = combineHapSeq(caseSeq);
    refGenoSeq = combineHapSeq(refSeq);
    testGenoSeq = combineHapSeq(refSeq);
    
    caseR = corrcoef(caseSeq);
    caseR(isnan(caseR)) = 0;
    refR = corrcoef(refSeq);
    refR(isnan(refR)) = 0;

    %calculate single allele frequence
    caseP = sum(caseSeq)/caseSize/2;
    refP = sum(refSeq)/refSize/2;
    testP = sum(testSeq)/testSize/2;
    
    %find sign rate and power
    
    %calculate Homer's attack's identification rate
    Tp_case = getTpM(caseGenoSeq, caseP, refP);
    Tp_test = getTpM(testGenoSeq, caseP, refP);
    zp = getThreshold(Tp_test, FDR);%get percentile
    idrP(i) = sum(Tp_case > zp);%get identification rate
    
    [idrSignPower(:,i) T(:,i)] = signRateIdr(caseSeq, testSeq, FDR, levels, caseR, refR);
        
    %calculate Yong's attack using copy sign from reference calculate R
    caseRCopySign = abs(caseR).*sign(refR);
    idrCS(i) = getIdr(caseSeq, testSeq, FDR, caseRCopySign, refR);
    
    if useEstR == 1
        %calculate estimate R
        caseEstR = estimateR(caseGenoSeq);
        caseEstR(isnan(caseEstR)) = 0;
        refEstR = estimateR(refGenoSeq);
        refEstR(isnan(refEstR)) = 0;
        idrEstSignPower(:,i) = signRateIdr(caseSeq, testSeq, FDR, levels, caseEstR, refEstR);
        caseEstRCopySign = abs(caseEstR).*sign(refR);
        idrCSEst(i) = getIdr(caseSeq, testSeq, FDR, caseEstRCopySign, refR);
    end
end
%profile viewer;

%% save
% configStr = ['case', num2str(caseSize),'ref', num2str(refSize), ...
%     'test', num2str(testSize), 'fdr', num2str(FDR), 'trial', num2str(trial), 'nSnps', num2str(nSnps), 'EstR', num2str(useEstR)];

fileName = [configStr, 'sp.mat'];
fprintf(1, 'write to %s\n', fileName);
save(fileName);

%% print result
fprintf(1, '%s\n', configStr);
fprintf(1, '\tHomer: \t%f\n', mean(idrP));
fprintf(1, '\trealR cp sign: \t%f\n', mean(idrCS));
fprintf(1, '\testR cp sign: \t%f\n', mean(idrCSEst));
fprintf(1, '\tmax Est R: \t%f\n', max(mean(idrEstSignPower,2)));
fprintf(1, '\tmax Real R: \t%f\n', max(mean(idrSignPower, 2)));

%% plot
h = figure;
maxT = 1;
line([0, maxT], [mean(idrP), mean(idrP)], 'Color', 'red', 'Marker', '.');
hold on;
line([0, maxT], [mean(idrCS), mean(idrCS)], 'Color', 'yellow', 'Marker', 'x');
line([0, maxT], [mean(idrCSEst), mean(idrCSEst)], 'Color', 'green', 'Marker', '.');
x = 1/levels:1/levels:1;
%x = mean(T, 2);
%set(gca,'XDir','reverse')
plot(x, mean(idrSignPower, 2), 'go-');
plot(x, mean(idrEstSignPower, 2), 'bx-');
legend('Homer', 'calcR cp refCalcR sign', 'calcEstR refCalcR sign','max power', 'idrEstSignPower', 2);
title(configStr);
ylabel('identification rate');
xlabel('level')
mkdir('./signPower');
saveas(h, ['./signPower/', configStr, '.pdf']);
hold off;
hasT = exist('T', 'var') > 0;
if hasT
    figure;
    h2 = plot(x, mean(T, 2));
    saveas(h2, ['./signPower/', configStr, 'T.pdf']);
end


