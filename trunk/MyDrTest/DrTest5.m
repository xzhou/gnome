function [StatS, StatR, StatT, Truth] = DrTest5()
% use outgroup as reference
% cd 'E:\1.2.2.Research Proteomics\GWAS information security\data\Real Dataset 2\Solver results from Rui\';
cd 'E:\1.2.2.Research Proteomics\GWAS information security\data\Hapmap chr10\solver result from Rui';
precision = 4;
Trials = 10;
estimateStd = 0;
Trialsonce = 1000;
if estimateStd == 0
    Trials = Trialsonce;
end
mini = 0;
adjustbias = 0 ; % estimate and adjust the bias from reference
%% read the haplotype sequences
% cd E:\1.2.2.Research Proteomics\GWAS information
% security\Mywork\five_population
%fastafile = 'chr10_FGFR2_200kb_phased_CEU.fasta';
fastafile_case = 'case.fasta';
fastafile_ctrl = 'ctrl.fasta';
fastafile_test = 'evolver.fasta';
% fastafile_test = 'init.fasta';
groupname = 'yri';
% outgroupname = 'yri';
outgroupname = 'jpt+chb';

% fastafile = 'chr10_FGFR2_200kb_phased_yri.fasta';
% fastafile_testgroup = 'chr10_FGFR2_200kb_phased_jpt+chb.fasta';
% groupname = 'yri';
% outgroupname = 'jpt+chb';

seqS = fastaread(fastafile_case);
seqR = fastaread(fastafile_ctrl); % test group
seqT = fastaread(fastafile_test);
nS = length(seqS);
nR = length(seqR);
nT = length(seqT);
Len = length(seqS(1).Sequence);
Len = 80;
%% divide sample into 2 groups: Sample group / reference group
int4S = zeros(nS,Len);
int4R = zeros(nR,Len);
int4T = zeros(nT,Len); % test group

for i = 1:nS
    seqS(i).Sequence = seqS(i).Sequence(1:Len);
    int4S(i,:)=nt2int(seqS(i).Sequence) - 1;
end
for i = 1:nR
    seqR(i).Sequence = seqR(i).Sequence(1:Len);
    int4R(i,:)=nt2int(seqR(i).Sequence) - 1;
end
for i = 1:nT
    seqT(i).Sequence = seqT(i).Sequence(1:Len);
    int4T(i,:)=nt2int(seqT(i).Sequence) - 1;
end

%% convert SNP nucleotide to 0/1
% define allele as allele 1 for each of the SNP location, thus that SNP
% with non 0/1 label can be converted to 0/1 label
% int4* is the 4 integer labeling of SNPs ,and int2* is the 2 integer
% labeling of SNPs
allele1 = int4T(end,:);
int2S = (int4S == repmat(allele1,nS,1)) + 0; 
int2R = (int4R == repmat(allele1,nR,1)) + 0;
int2T = (int4T == repmat(allele1,nT,1)) + 0;


%% compute the r-values of the sample group, and reference group
all_r_S = corrcoef(int2S);
all_r_S(isnan(all_r_S))= 0; % deal with invariant sites
all_r_R = corrcoef(int2R);
all_r_R(isnan(all_r_R))= 0; % deal with invariant sites

%% compute the Tr test statistic values
%% simulation using reference group and obtain the estimated variance
%% using z-test to optain the p-values for each sample/reference individual
StatS.Tr = zeros(nS,1);
StatS.std = StatS.Tr;
StatS.p = StatS.Tr;
StatR.Tr = zeros(nR,1);
StatR.std = StatR.Tr;
StatR.p = StatR.Tr;
StatT.Tr = zeros(nT,1);
StatT.std = StatT.Tr;
StatT.p = StatT.Tr;

iMCmodel_R = iMCmodelBuild(seqR, 0); % use refernence to build a model for simulation
Truth = zeros(nS + nR + nT,1);
std = getstdTr([], iMCmodel_R, allele1, nS, nR, Trialsonce, precision)
for i = 1:nS
    StatS.Tr(i) = getTr(int2S(i,:), all_r_S, all_r_R, adjustbias);
    if estimateStd ==1 
        StatS.std(i) = getstdTr(int4S(i,:), iMCmodel_R, allele1, nS, nR, Trials, precision);
    end
    sim = sum(repmat(int4S(i,:),nS,1)==int4S,2);
    Truthmax(i) = 1;
    Truthmean(i) = mean(sim)/Len;
end

for i = 1:nR
    StatR.Tr(i) = getTr(int2R(i,:), all_r_S, all_r_R, adjustbias);
    if estimateStd ==1 
        StatR.std(i) = getstdTr(int4R(i,:), iMCmodel_R, allele1, nS, nR, Trials, precision);
    end
    sim = sum(repmat(int4R(i,:),nS,1)==int4S,2);
    Truthmax(i + nS) = max(sim)/Len;
    Truthmean(i + nS) = mean(sim)/Len;
end

for i = 1:nT
    StatT.Tr(i) = getTr(int2T(i,:), all_r_S, all_r_R, adjustbias);
    if estimateStd ==1 
        StatT.std(i) = getstdTr(int4T(i,:), iMCmodel_R, allele1, nS, nR, Trials, precision);
    end
    sim = sum(repmat(int4T(i,:),nS,1)==int4S,2);
    Truthmax(i + nS + nR) = max(sim)/Len;
    Truthmean(i + nS + nR) = mean(sim)/Len;
end
StatS.Tr = StatS.Tr/sqrt(Len*(Len-1)/2);
StatR.Tr = StatR.Tr/sqrt(Len*(Len-1)/2);
StatT.Tr = StatT.Tr/sqrt(Len*(Len-1)/2);
if estimateStd == 0 
    StatS.std = StatS.std+std;
    StatR.std = StatR.std+std;
    StatT.std = StatT.std+std;
end
StatS.p = normcdf(-StatS.Tr,0,StatS.std);
StatR.p = normcdf(-StatR.Tr,0,StatR.std);
StatT.p = normcdf(-StatT.Tr,0,StatT.std);
% plot the p-values

indxSim95 = find(Truthmax>0.95 & Truthmax<1); 
indxSim95 = indxSim95(find(indxSim95>nS));
indxSim = find(Truthmax==1); 
indxSim = indxSim(find(indxSim>nS));
Trall = [StatS.Tr; StatR.Tr; StatT.Tr];
Stdall = [StatS.std; StatR.std; StatT.std];
pall = [StatS.p; StatR.p; StatT.p];
indxS = 1:nS;
indxR = nS+1:nS+nR;
indxT = nS+nR+1:nS+nR+nT;
indxR = setdiff(indxR,union(indxSim95, indxSim));
indxT = setdiff(indxT,union(indxSim95, indxSim));

figure;
hold on;
title(['Tr values; ' 'precision: ' num2str(precision) '; Trials: ' num2str(Trials)]);
plot(indxS, Trall(indxS), '.r');
plot(indxR, Trall(indxR), '.g');
plot(indxT,  Trall(indxT), '.b');
plot(indxSim95, Trall(indxSim95), 'xr');
plot(indxSim, Trall(indxSim), '*r');
xlabel('indx of individuals');
ylabel('Tr values');
legend({'Sample individuals' 'Reference individuals' ...
    'solver results' 'Simmax>0.95'...
    'Simmax==1'});

figure;
hold on;
title(['Std values; ' 'precision: ' num2str(precision) '; Trials: ' num2str(Trials)]);
plot(indxS, Stdall(indxS), '.r');
plot(indxR, Stdall(indxR), '.g');
plot(indxT,  Stdall(indxT), '.b');
plot(indxSim95, Stdall(indxSim95), 'xr');
plot(indxSim, Stdall(indxSim), '*r');
xlabel('indx of individuals');
ylabel('Std values');
legend({'Sample individuals' 'Reference individuals' ...
    'solver results' 'Simmax>0.95'...
    'Simmax==1'});

hold off;
figure;
semilogy(indxS, pall(indxS), '.r');
hold on;
semilogy(indxR, pall(indxR), '.g');
semilogy(indxT, pall(indxT), '.b');
semilogy(indxSim95, pall(indxSim95), 'xr');
semilogy(indxSim, pall(indxSim), '*r');
semilogy([1 nS+nR+nT], [1e-2 1e-2], '--b');
semilogy([1 nS+nR+nT], [1e-5 1e-5], '--r');
xlabel('indx of individuals');
ylabel('p value')
title(['p values; ' 'precision: ' num2str(precision) '; Trials: ' num2str(Trials)]);
legend({'Sample individuals' 'Reference individuals' ...
    'solver results' 'Simmax>0.95'...
    'Simmax==1' 'p-value 0.01 cutoff' 'p-value 1.0E-5 cutoff'});

figure;
hold off;
semilogy(Truthmax(indxS), pall(indxS), '.r');
hold on;
semilogy(Truthmax(indxR), pall(indxR), '.g');
semilogy(Truthmax(indxT), pall(indxT), '.b');
xlabel('indx of individuals');
ylabel('p value');
xlabel('max similarity to the sample individuals');
title(['max similarity vs. p-value values; ' 'precision: ' num2str(precision) '; Trials: ' num2str(Trials)]);
legend({'Sample individuals' 'Reference individuals' ...
    'solver results'});

figure;
hold off;
semilogy(Truthmean(indxS), pall(indxS), '.r');
hold on;
semilogy(Truthmean(indxR), pall(indxR), '.g');
semilogy(Truthmean(indxT), pall(indxT), '.b');
ylabel('p value');
xlabel('mean similarity to the sample individuals');
title(['mean similarity vs. p-value values; ' 'precision: ' num2str(precision) '; Trials: ' num2str(Trials)]);
legend({'Sample individuals' 'Reference individuals' 'solver results'});

nIDS_e2 = sum(StatS.p<0.01)
nIDR_e2 = sum(StatR.p<0.01)
nIDT_e2 =sum(StatT.p<0.01)

nIDS_e10 = sum(StatS.p<10e-10)
nIDR_e10 = sum(StatR.p<10e-10)

nIDT_e10 = sum(StatT.p<10e-10)

function [std_Dr, mean_Dr] = getstdTr(int4, iMCmodel, refhaplotype, nS, nR, Trials, precision)
Len = length(int4);
if Len ~= 0
    allele1 = int4;
else
    allele1 = iMCgenerate(iMCmodel,1);
end
Len = length(allele1);
% Sample one individual form the model
A = (allele1 == refhaplotype) + 0;
A2 = (2*A'-1)*(2*A-1);

D_r_M0 = zeros(Trials,1);

for i = 1: Trials
    sample_R = iMCgenerate(iMCmodel,nR);
    sample_S = iMCgenerate(iMCmodel,nS);
    sample_R = (sample_R == repmat(refhaplotype,nR,1)) + 0;
    sample_S = (sample_S == repmat(refhaplotype,nS,1)) + 0;

%     sample_S(end,:)=A;
    all_r_R = corrcoef(sample_R);
    all_r_R(isnan(all_r_R))= 0; % deal with invariant sites
    all_r_S = corrcoef(sample_S);
    all_r_S(isnan(all_r_S))= 0; % deal with invariant sites

    % Test statistic for r
    D_r_S(i) = sum(sum((all_r_S - all_r_R).* A2))/2;
end
D_r_S = D_r_S/sqrt(Len*(Len-1)/2);
figure
hist(D_r_S);
mean_Dr = mean(D_r_S);
std_Dr = sqrt(var(D_r_S));
return

function Tr = getTr(Y, r_S, r_R, adjustbias)
A2 = (2*Y'-1)*(2*Y-1);
if adjustbias == 1
    A2 = A2 - 1;
end
Tr = sum(sum((r_S - r_R).* A2))/2;

return