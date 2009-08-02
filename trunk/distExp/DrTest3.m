function [StatS, StatR, StatT, Truth] = DrTest3()
% do the test given r^2 from snp.plotter for the target sample and hapmap
% reference
precision = 4;
Trials = 10;
estimateStd = 0;
Trialsonce = 1000;
if estimateStd == 0
    Trials = Trialsonce;
end
signByRef = 0; % use reference to estimate sign
mini = 0;
fromSnpPlotter = 1; % use snp plotter r values
percentage_keep = 1; % only keep the x percentage of the lower values
withmissingSNP = 1; % decided the missing SNPs

%% read the haplotype sequences
% cd E:\1.2.2.Research Proteomics\GWAS information security\Mywork\five_population
%% 
% r values from snp.plotter
%% To change: input files
% cd 'E:\1.2.2.Research_Proteomics\GWAS information security\data\Real Dataset 2\Statistical test'

Sample_snpplotterFile = '80SNP_CEU_sim_4000seq_sample1000.ped.freq.txt';
Sample_haplotypeFile = '80SNP_CEU_sim_4000seq_sample1000.fasta';

Ref_snpplotterFile = '80SNP_CEU_sim_4000seq_control1000.freq.txt';
Ref_haplotypeFile = '80SNP_CEU_sim_4000seq_control1000.fasta';

Test_haplotypeFile = 'hapmap_chr7_80SNP_CEU_haplotype_83_164.fasta';
Out_haplotypeFile = 'hapmap_chr7_80SNP_YRI_haplotype.fasta';

%%
[R2 R1]= readSNPplotter(Ref_snpplotterFile);
[S2 S1]= readSNPplotter(Sample_snpplotterFile);
LenS = size(S1.p,2);
LenR = size(R1.p,2);
% missing snps in hapmap phasing
Smissingcol = [];
Phasingmiss = [];
Plottermiss = [];
if withmissingSNP == 1
    SNPidGenometype = getID('SNPID80.txt');
    SNPidHapmapPhased = getID('SNPfound.list');
    [tmp Phasingmiss] = setdiff(SNPidGenometype,SNPidHapmapPhased);

    % misssing snps in snp.plotter output with hapmap genotype as input
    Plottermiss = [17 46]; % these SNPs are missed in the reference data
    % Plottermiss =  [17 42 46];
    % handle the missing cols
    Smissingcol = union(Phasingmiss, Plottermiss);
end

Sleftcol = setdiff(1:LenS, Smissingcol);
Len = length(Sleftcol); % final length
S2clean.r(1:Len,1:Len) =  S2.r(Sleftcol,Sleftcol);
S1clean.name =  {S1.name{1,Sleftcol}; S1.name{2,Sleftcol}};
S1clean.p =  S1.p(Sleftcol);

Rmissingcol = indxshift(Phasingmiss, Plottermiss);
leftcol = setdiff(1:LenR, Rmissingcol);
Len = length(leftcol); % final length
R2clean.r(1:Len,1:Len) =  R2.r(leftcol,leftcol);
R1clean.name =  {R1.name{1,leftcol}; R1.name{2,leftcol}};
R1clean.p =  R1.p(leftcol);

[S2clean.r, S1clean] = alleleMatch12(S2clean.r, S1clean, R2clean.r, R1clean);



SeqS = fastaread(Sample_haplotypeFile);
SeqR = fastaread(Ref_haplotypeFile);
SeqOther = fastaread(Test_haplotypeFile); %ingroup test
SeqT = fastaread(Out_haplotypeFile); % outgroup as test 
groupname = 'ceu';
outgroupname = 'yri';
SeqT = [SeqOther; SeqT];  % test group

%% reference haplotype sequences
Rmissingcol = indxshift(Plottermiss, Phasingmiss);
LenR = length(SeqR(1).Sequence);
leftcol = setdiff(1:LenR, Rmissingcol);
Len = length(leftcol);
nS = length(SeqS);
nR = length(SeqR);
nT = length(SeqT);
int4S = zeros(nS,Len);
int4R = zeros(nR,Len);
int4T = zeros(nT,Len); % test group
for i = 1:length(SeqR)
    SeqR(i).Sequence = SeqR(i).Sequence(leftcol);
    int4R(i,:) = nt2int(SeqR(i).Sequence) - 1;
end
for i = 1:nT   
    SeqT(i).Sequence = SeqT(i).Sequence(leftcol);
    int4T(i,:) = nt2int(SeqT(i).Sequence) - 1;
end
for i = 1:length(SeqS)
    i;
    seq = SeqS(i).Sequence(Sleftcol);
    SeqS(i).Sequence = alleleTranslate(seq, S1clean.oldname, S1clean.name);
    int4S(i,:) = nt2int(SeqS(i).Sequence)-1;
end


%% convert SNP nucleotide to 0/1
% define allele as allele 1 for each of the SNP location, thus that SNP
% with non 0/1 label can be converted to 0/1 label
% int4* is the 4 integer labeling of SNPs ,and int2* is the 2 integer
% labeling of SNPs
seq = repmat('1',1,Len);
seq = alleleTranslate(seq, S1clean.oldname, S1clean.name)
allele1 = nt2int(seq)-1;
int2S = (int4S == repmat(allele1,nS,1)) + 0; 
int2R = (int4R == repmat(allele1,nR,1)) + 0;
int2T = (int4T == repmat(allele1,nT,1)) + 0;

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

iMCmodel_R = iMCmodelBuild(SeqR, 0); % use refernence to build a model for simulation
Truth = zeros(nS + nR + nT,1);
std = getstdTr([], iMCmodel_R, allele1, nS, nR, Trialsonce, precision)
all_r_S = S2clean.r;
all_r_R = R2clean.r;

%% compute the r-values of the sample group, and reference group
if fromSnpPlotter == 0
    all_r_S = corrcoef(int2S);
    all_r_S(isnan(all_r_S))= 0; % deal with invariant sites
    all_r_R = corrcoef(int2R);
    all_r_R(isnan(all_r_R))= 0; % deal with invariant sites
end

(sum(sum(sign(all_r_S)== sign(all_r_R))))/(prod(size(all_r_R)))
x = reshape(all_r_R,1, prod(size(all_r_R)));
x = sort(abs(x));
indx = floor(length(x)*percentage_keep);
threshold = x(indx)
indxall = (abs(all_r_S))> threshold;
all_r_S(indxall) = 0;
all_r_R(indxall) = 0;    

if signByRef == 1
    all_r_S = abs(all_r_S).*sign(all_r_R);
end
% indx = (sign(all_r_S)~=sign(all_r_R));
% all_r_S(indx)=0;
% all_r_R(indx)=0;

for i = 1:nS
    StatS.Tr(i) = getTr(int2S(i,:), all_r_S, all_r_R);
    % StatR.Tr(i) = getTr(int2R(i,:), all_r_S, all_r_R);
    if estimateStd ==1 
        StatS.std(i) = getstdTr(int4S(i,:), iMCmodel_R, allele1, nS, nR, Trials, precision);
        % StatR.std(i) = getstdTr(int4R(i,:), iMCmodel_R, allele1, nS, nR, Trials, precision);
    end
    sim = sum(repmat(int4S(i,:),nS,1)==int4S,2);
    Truthmax(i) = 1;
    Truthmean(i) = mean(sim)/Len;
%     sim = sum(repmat(int4R(i,:),nS,1)==int4S,2);
%     Truthmax(i + nS) = max(sim)/Len;
%     Truthmean(i + nS) = mean(sim)/Len;
end

for i = 1:nR
    % StatS.Tr(i) = getTr(int2S(i,:), all_r_S, all_r_R);
    StatR.Tr(i) = getTr(int2R(i,:), all_r_S, all_r_R);
    if estimateStd ==1 
        % StatS.std(i) = getstdTr(int4S(i,:), iMCmodel_R, allele1, nS, nR, Trials, precision);
        StatR.std(i) = getstdTr(int4R(i,:), iMCmodel_R, allele1, nS, nR, Trials, precision);
    end
%     sim = sum(repmat(int4S(i,:),nS,1)==int4S,2);
%     Truthmax(i) = 1;
%     Truthmean(i) = mean(sim)/Len;
    sim = sum(repmat(int4R(i,:),nS,1)==int4S,2);
    Truthmax(i + nS) = max(sim)/Len;
    Truthmean(i + nS) = mean(sim)/Len;
end

for i = 1:nT
    StatT.Tr(i) = getTr(int2T(i,:), all_r_S, all_r_R);
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
nSample = length(SeqS)+length(SeqR)+length(SeqOther);
indxS = 1:nS;
indxR = nS+1:nS+nR;
indxOther = nS+nR+1:nSample;
indxT = nSample+1:nS+nR+nT;
indxR = setdiff(indxR,union(indxSim95, indxSim));
indxOther = setdiff(indxOther,union(indxSim95, indxSim));
indxT = setdiff(indxT,union(indxSim95, indxSim));
figure;
hold on;
title(['Tr values; ' 'precision: ' num2str(precision) '; Trials: ' num2str(Trials)]);
plot(indxS, Trall(indxS), '.r');
plot(indxR, Trall(indxR), '.g');
plot(indxOther,  Trall(indxOther), '.b');
plot(indxT,  Trall(indxT), '.k');
plot(indxSim95, Trall(indxSim95), 'xr');
plot(indxSim, Trall(indxSim), '*r');
xlabel('indx of individuals');
ylabel('Tr values');
legend({[groupname ' Sample individuals'] [groupname ' Reference individuals']...
    [groupname ' Other individuals'] [outgroupname ' Out group']  'Simmax>0.95'...
    'Simmax==1'});
if fromSnpPlotter == 0
    saveas(gcf,'Tr_80SNP_dataset_yri_as_outgroup(r values from Haplotype).fig'); 
    saveas(gcf,'Tr_80SNP_dataset_yri_as_outgroup(r values from Haplotype).png');
else 
    saveas(gcf,'Tr_values_80SNP_dataset_yri_as_outgroup (r values from snp.plotter).fig'); 
    saveas(gcf,'Tr_values_80SNP_dataset_yri_as_outgroup (r values from snp.plotter).png');
end

figure;
hold on;
title(['Std values; ' 'precision: ' num2str(precision) '; Trials: ' num2str(Trials)]);
plot(indxS, Stdall(indxS), '.r');
plot(indxR, Stdall(indxR), '.g');
plot(indxOther,  Stdall(indxOther), '.b');
plot(indxT,  Stdall(indxT), '.k');
plot(indxSim95, Stdall(indxSim95), 'xr');
plot(indxSim, Stdall(indxSim), '*r');
xlabel('indx of individuals');
ylabel('Std values');
legend({[groupname ' Sample individuals'] [groupname ' Reference individuals']...
    [groupname ' Other individuals'] [outgroupname ' Out group']  'Simmax>0.95'...
    'Simmax==1'});

hold off;
figure;
semilogy(indxS, pall(indxS), '.r');
hold on;
semilogy(indxR, pall(indxR), '.g');
semilogy(indxOther, pall(indxOther), '.b');
semilogy(indxT,  pall(indxT), '.k');
semilogy(indxSim95, pall(indxSim95), 'xr');
semilogy(indxSim, pall(indxSim), '*r');
semilogy([1 nS+nR+nT], [1e-2 1e-2], '--b');
semilogy([1 nS+nR+nT], [1e-5 1e-5], '--r');
xlabel('indx of individuals');
ylabel('p value')
lgd = {[groupname ' Sample individuals'] [groupname ' Reference individuals']...
    [groupname ' Other individuals'] [outgroupname ' Out group']};
if length(indxSim95)>0
    lgd{end+1} = 'Simmax>0.95';
end
if length(indxSim)>0
    lgd{end+1} = 'Simmax==1';
end
lgd = {lgd{1:end} 'p-value 0.01 cutoff' 'p-value 1.0E-5 cutoff'};
title(['p values; ' 'precision: ' num2str(precision) '; Trials: ' num2str(Trials)]);
legend(lgd);
if fromSnpPlotter == 0
    saveas(gcf,'p_values_80SNP_dataset_yri_as_outgroup (r values from Haplotype).fig'); 
    saveas(gcf,'p_values_80SNP_dataset_yri_as_outgroup (r values from Haplotype).png');
else 
    saveas(gcf,'p_values_80SNP_dataset_yri_as_outgroup (r values from snp.plotter).fig'); 
    saveas(gcf,'p_values_80SNP_dataset_yri_as_outgroup (r values from snp.plotter).png');
end

figure;
hold off;
semilogy(Truthmax(1:nS), StatS.p, '.r');
hold on;
semilogy(Truthmax(nS+1:nS+nR), StatR.p, '.g');
semilogy(Truthmax(nS+nR+1:nSample), StatT.p(1:nSample-(nS+nR)), '.b');
semilogy(Truthmax(nSample+1:nS+nR+nT), StatT.p(nSample-(nS+nR)+1:nT), '.k');
xlabel('indx of individuals');
ylabel('p value');
xlabel('max similarity to the sample individuals');
title(['max similarity vs. p-value values; ' 'precision: ' num2str(precision) '; Trials: ' num2str(Trials)]);
legend({[groupname ' Sample individuals'] [groupname ' Reference individuals'] [groupname ' Other individuals'] [outgroupname ' Out group']  ['Simmax>0.95']});

figure;
hold off;
semilogy(Truthmean(1:nS), StatS.p, '.r');
hold on;
semilogy(Truthmean(nS+1:nS+nR), StatR.p, '.g');
semilogy(Truthmean(nS+nR+1:nSample), StatT.p(1:nSample-(nS+nR)), '.b');
semilogy(Truthmean(nSample+1:nS+nR+nT), StatT.p(nSample-(nS+nR)+1:nT), '.k');
ylabel('p value');
xlabel('mean similarity to the sample individuals');
title(['mean similarity vs. p-value values; ' 'precision: ' num2str(precision) '; Trials: ' num2str(Trials)]);
legend({[groupname ' Sample individuals'] [groupname ' Reference individuals'] [groupname ' Other individuals'] [outgroupname ' Out group']});

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

function Tr = getTr(Y, r_S, r_R)
A2 = (2*Y'-1)*(2*Y-1);
Tr = sum(sum((r_S - r_R).* A2))/2;
return
