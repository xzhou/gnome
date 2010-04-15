function out = DrTest(fastafile, threshold, precision, N)
startParallel();
cd '/home/xzhou/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/fastPhase';
cd '/home/xzhou/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/fastPhase523';
nP = 100;
nM = 100;
if nargin >=4
    nP = N;
    nM = N;
end
Trials = 1000;
Len = 100;
flag_usereal = 1;
degFreedomFlag = 0; % estimate the degree of freedom of the haplotypes
fastafile = 'Affx_gt_58C_Chiamo_07.tped.200snp.extract.inp.fasta';
fastafile = 'Affx_gt_58C_Chiamo_07.tped.600SNP.extract.inp.fasta';
% fastafile = 'chr10_FGFR2_200kb_phased_yri.fasta';
frompairwise = 1; % compute the statistics from pairwise freq
% Tr_Tp = 1; % use Tr+Tp instead of Tr
estimateSign = 0; % estimate the sign of r in the target sample from the reference
removebySign = 0; % 1, remove the r values with the same sign in sample and reference.
% -1, remove the r values with different sign in sample
% and references
% else, do not remove r values by sign
partialReverseSign = 0; % partially reverse the sign of the "removebySign" part; or totaly remove the "removebySign" part
pseudocount = 1;
if nargin < 2 || length(threshold)==0
    threshold = 0.1; %0.5, 0.2, 0.1, 0.01, 0.001
end
removeSmallr = 0;
flag_twoPopulation = 0; % two models will be simulated for the references and the sample
if nargin < 3 || length(precision)==0
    precision = 10; %1, 2, 3, 4, 5
end

%% read data from fasta file
rawFastaData = fastaread(fastafile);
[hap01Seq, alleleMapping] = encodeRawFastaSeq(rawFastaData);
[nInvariantSnps, indexs] = findInvariantSnps(hap01Seq);
hap01Seq(:, indexs) = [];%remove invariant snps
alleleMapping(:, indexs) = [];
hap01Seq = unique(hap01Seq, 'rows');
[nUnique realLen] = size(hap01Seq);
if realLen > Len
    hap01Seq = hap01Seq(:, 1:Len);
else
    warning('x:y', 'not enough snp, using all available');
    Len = realLen;
end

allele1 = hap01Seq(1,:);

if degFreedomFlag ==1
    d = zeros(200,1);
    for i = 1:200
        d(i) = sum(iMCgenerate(iMCmodel,1)==iMCgenerate(iMCmodel,1));
    end
    mean_d = mean(d);
    effictive_deg_freedom = (Len - mean_d)*2;
end

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock))); 

% Sample one individual form the model,
[A, idx] = randomKRow(hap01Seq, 1);
hap01Seq(idx, :) = [];%remove A;
A = (A == allele1) + 0;
A1 = (2*A'-1);
A2 = (2*A'-1)*(2*A-1);  %sign

% this is the test statistics each round
D_r_M0 = zeros(Trials,1);
D_r_MA = zeros(Trials,1);
D_p_M0 = zeros(Trials,1);
D_p_MA = zeros(Trials,1);

countsmall = 0;
if flag_twoPopulation ==1
    iMCmodel1 = iMCmodelBuild(iMCgenerate(iMCmodel,50),0);
    iMCmodel2 = iMCmodelBuild(iMCgenerate(iMCmodel,50),0);
end

%% allocate space for test values
allTcMA = zeros(Len, Len, Trials);
allTcM0 = zeros(Len, Len, Trials);

%% repeat for 10000 times and calculate the distribution of r_M0 r_MA
%profile on;
for i = 1:Trials
    fprintf(1, '%d\n', i);
    
    hap01SeqTemp = hap01Seq;%make a copy
    [sample_P, idx] = randomKRow(hap01SeqTemp, nP);
    hap01SeqTemp(idx, :) = [];
    [sample_M, ~] = randomKRow(hap01SeqTemp, nM);

    %encode
    sample_P = (sample_P == repmat(allele1,nP,1)) + 0;
    sample_M = (sample_M == repmat(allele1,nM,1)) + 0;
    
    all_r_P = corrcoef(sample_P);
    all_r_P(isnan(all_r_P))= 0; % deal with invariant sites
    
    all_r_M0 = corrcoef(sample_M);
    all_r_M0(isnan(all_r_M0))= 0; % deal with invariant sites
    
    all_p_P = sum(sample_P,1)/nP;
    all_p_M0 = sum(sample_M,1)/nM;
    
    if frompairwise == 1
        [Tc_M0(i) Tc_log_M0(i) Tr_M0(i) TD_M0(i) allTcM0(:, :, i)] = TcRetAll(sample_M, sample_P, A);
        sample_M(end,:)=A;
        [Tc_MA(i) Tc_log_MA(i) Tr_MA(i) TD_MA(i) allTcMA(:, :, i)] = TcRetAll(sample_M, sample_P, A);
    else
        sample_M(end,:)=A;
    end
    all_r_MA = corrcoef(sample_M);
    all_r_MA(isnan(all_r_MA))= 0; % deal with invariant sites
    all_p_MA = sum(sample_M,1)/nM;
    
    all_r_P = round(all_r_P*10^precision)./10^precision;
    all_r_M0 = round(all_r_M0*10^precision)./10^precision;
    all_r_MA = round(all_r_MA*10^precision)./10^precision;
    all_p_P = round(all_p_P*10^precision)./10^precision;
    all_p_M0 = round(all_p_M0*10^precision)./10^precision;
    all_p_MA = round(all_p_MA*10^precision)./10^precision;
    
    % sign_agree_MAM0 = sum(sum(sign(all_r_MA)==sign(all_r_M0)))/prod(size(all_r_M0))
    % Test statistic for r
    % just search for the flag "estimateSign" you will see the relavent
    % block of code --- the signs here were repalced by the signs of r computed from the reference sample.
    if estimateSign == 1
        sign_agree_M0P = sum(sum(sign(all_r_M0)==sign(all_r_P)))/prod(size(all_r_M0));
        sign_agree_MAP = sum(sum(sign(all_r_MA)==sign(all_r_P)))/prod(size(all_r_M0));
        sign_agree_MAM0 = sum(sum(sign(all_r_MA)==sign(all_r_M0)))/prod(size(all_r_M0));
        
        %replace sign here
        signMatrix = sign(all_r_P);
        
        %all_r_M0 = abs(all_r_M0).*sign(all_r_P);
        %all_r_MA = abs(all_r_MA).*sign(all_r_P);
        
        %all_r_M0 = abs(all_r_M0).*signMatrix;
        %all_r_MA = abs(all_r_MA).*signMatrix;
        
        %modified by Xiaoyong
        all_r_M0 = replaceSignByT(all_r_M0, 0.1);
        all_r_MA = replaceSignByT(all_r_MA, 0.1);
    end
    
    
    indx0 = [];
    indxA = [];
    if removebySign == 1
        indx0 = sign(all_r_M0)==sign(all_r_P);
        indxA = sign(all_r_MA)==sign(all_r_P);
    elseif removebySign == -1
        indx0 = sign(all_r_M0)==sign(all_r_P);
        indxA = sign(all_r_MA)==sign(all_r_P);
    end
    if partialReverseSign == 1
        all_r_M0(indx0) = -all_r_M0(indx0);
        all_r_MA(indxA) = -all_r_MA(indxA);
    else
        all_r_M0(indx0) = all_r_P(indx0);
        all_r_MA(indxA) = all_r_P(indxA);
    end
    
    
    if removeSmallr == 1
        countsmall = countsmall + sum(sum(abs(all_r_M0)<threshold));
        %         r = all_r_M0(all_r_M0>=0);
        %         logr = log10(r);
        %         logr(~isfinite(logr))=1;
        %         hist(logr)
        all_r_M0(abs(all_r_M0)<threshold) = 0;
        all_r_MA(abs(all_r_MA)<threshold) = 0;
        all_r_P(abs(all_r_P)<threshold) = 0;
    end
    
    %calculate the pairwise test
    D_r_M0(i) = sum(sum((all_r_M0 - all_r_P).* A2))/2;
    D_r_MA(i) = sum(sum((all_r_MA - all_r_P).* A2))/2;
    
    % adjust for systemetic bias
    % D_radjust_M0(i) = D_r_M0(i) - sum(sum(all_r_M0 - all_r_P))/2;
    % D_radjust_MA(i) = D_r_MA(i) - sum(sum(all_r_MA - all_r_P))/2;
    
    %     D_radjust2_M0(i) = D_r_M0(i)- sum(sum(all_r_M0 - all_r_P))/2;
    %     D_radjust2_MA(i) = D_r_MA(i) - sum(sum(all_r_MA - all_r_P))/2;
    %
    %     D_rsign_M0(i) = sum(sum(sign(all_r_M0 - all_r_P).* A2))/2;
    %     D_rsign_MA(i) = sum(sum(sign(all_r_MA - all_r_P).* A2))/2;
    
    % correct for the variance of (all_r_M0 - all_r_P);
    % var(r|bivariate normal population) = (1-r^2)^2/(n-1)
    %     r_correct = (1-all_r_P.^2);
    %     r_correct(r_correct < 0.025)=1;
    %     D_rcorrection_M0(i) = sum(sum((all_r_M0 - all_r_P)./r_correct.* A2))/2;
    %     D_rcorrection_MA(i) = sum(sum((all_r_MA - all_r_P)./r_correct.* A2))/2;
    %
    % Test statistic for estimated frequency
    D_p_M0(i) = (all_p_M0 - all_p_P)* A1;
    D_p_MA(i) = (all_p_MA - all_p_P)* A1;
    %     D_psign_M0(i) = sign(all_p_M0 - all_p_P)* A1;
    %     D_psign_MA(i) = sign(all_p_MA - all_p_P)* A1;
    %     p_correct = (all_p_P.*(1-all_p_P));
    %     p_correct(p_correct<0.05)=1;
    %     D_pcorrection_M0(i) = (all_p_M0 - all_p_P)./p_correct* A1;
    %     D_pcorrection_MA(i) = (all_p_MA - all_p_P)./p_correct* A1;
end
save;
save(['real', num2str(Trials), 'nP', num2str(nP), '.mat']);


if frompairwise == 1
    Tc_MA = Tc_MA/sqrt(Len*(Len-1)/2);
    Tc_log_MA = Tc_log_MA/sqrt(Len*(Len-1)/2);
    Tc_M0 = Tc_M0/sqrt(Len*(Len-1)/2);
    Tr_MA = Tr_MA/sqrt(Len*(Len-1)/2);
    Tr_M0 = Tr_M0/sqrt(Len*(Len-1)/2);
    
    Tc_log_M0 = Tc_log_M0/sqrt(Len*(Len-1)/2);
    TD_MA = TD_MA/sqrt(Len*(Len-1)/2);
    TD_M0 = TD_M0/sqrt(Len*(Len-1)/2);
end

%% analysis of each pair
allTcM0 = diagZero(allTcM0);
allTcMA = diagZero(allTcMA);
meanTcM0 = mean(allTcM0, 3);
meanTcMA = mean(allTcMA, 3);
stdTcM0 = std(allTcM0, 3);
stdTcMA = std(allTcMA, 3);

uM0 = sum(meanTcM0);
uMA = sum(meanTcMA);

deltaM0 = sqrt(sum(stdTcM0.*stdTcM0,3));
deltaMA = sqrt(sum(stdTcMA.*stdTcMA,3));

%assuming standard deviation calculate the max power
q95 = norminv(0.95, uM0, deltaM0);
power95 = norminv(q95, uMA, deltaMA);
fprintf(1, 'MAX power %f', power95);


%%normalize the Tr test
D_r_M0 = D_r_M0/sqrt(Len*(Len-1)/2);
D_r_MA = D_r_MA/sqrt(Len*(Len-1)/2);
mean_r0 = mean(D_r_M0);
mean_rA = mean(D_r_MA);
var_r0 = var(D_r_M0);
var_rA = var(D_r_MA);

% D_radjust_M0 = D_radjust_M0/sqrt(Len*(Len-1)/2);
% D_radjust_MA = D_radjust_MA/sqrt(Len*(Len-1)/2);
% mean_radjust0 = mean(D_radjust_M0);
% mean_radjustA = mean(D_radjust_MA);
% var_radjust0 = var(D_radjust_M0);
% var_radjustA = var(D_radjust_MA);

%%normalize the Tp test
D_p_M0 = D_p_M0/sqrt(Len);
D_p_MA = D_p_MA/sqrt(Len);
% D_rsign_M0 = D_rsign_M0/sqrt(Len*(Len-1)/2);
% D_rsign_MA = D_rsign_MA/sqrt(Len*(Len-1)/2);
% D_psign_M0 = D_psign_M0/sqrt(Len);
% D_psign_MA = D_psign_MA/sqrt(Len);
% D_rcorrection_M0 = D_rcorrection_M0/sqrt(Len*(Len-1)/2);
% D_rcorrection_MA = D_rcorrection_MA/sqrt(Len*(Len-1)/2);
% mean_r0 = mean(D_rcorrection_M0);
% mean_rA = mean(D_rcorrection_MA);
% var_r0 = var(D_rcorrection_M0);
% var_rA = var(D_rcorrection_MA);
% D_pcorrection_M0 = D_pcorrection_M0/sqrt(Len);
% D_pcorrection_MA = D_pcorrection_MA/sqrt(Len);
countsmall = countsmall/Trials/((Len-1)*Len);

% close all;
% textTitle = {'Statistic for r, 2000 SNP, control size 100, model group size 50';'Null: M does not contain the individual'};

%% plot
textTitle = {'Tr test M0 MA'};
power_Dr=plot2hist(D_r_M0,D_r_MA,textTitle);
xlabel('T_r');
ylabel('Density of T_r');
title('')
fprintf(1, 'power Dr %f\n', power_Dr);

textTitle = {'Tr+Tp M0 MA'};
D_r_M0 = D_r_M0 + D_p_M0;
D_r_MA = D_r_MA + D_p_MA;
power_DrDp=plot2hist(D_r_M0,D_r_MA,textTitle);
xlabel('T_r + T_p');
ylabel('Density of T_r + T_p');
title('')

% textTitle = {['case:' num2str(nM) ' control:' num2str(nP) ' Tr, FGFR2 locus, 174 SNPs, control size 100, model group size 50'];'Null: M does not contain the individual'};
% power_Tr=plot2hist(Tr_M0,Tr_MA,textTitle);
% xlabel('T_r');
% ylabel('Density of T_r');
% title('')

% textTitle = {['case:' num2str(nM) ' control:' num2str(nP) ' TD, FGFR2 locus, 174 SNPs,'];'Null: M does not contain the individual'};
% power_TD = plot2hist(TD_M0,TD_MA,textTitle);

textTitle = {'Tp M0 MA'};
power_Dp = plot2hist(D_p_M0,D_p_MA,textTitle);
xlabel('T_p');
ylabel('Density of T_p');
title('')
fprintf(1, 'T_p = %f\n', power_Dp);

%     out = [threshold; countsmall; power_Dr; power_Dp; power_DrDp];
%     textTitle = {['case:' num2str(nM) ' control:' num2str(nP) ' Tc: ' fastafile];'Null: M does not contain the individual'};
%     plot2hist(Tc_M0,Tc_MA,textTitle)
%
textTitle = {['case:' num2str(nM) ' control:' num2str(nP) ' Tc\_log: ' fastafile];'Null: M does not contain the individual'};
power_Tc = plot2hist(Tc_log_M0,Tc_log_MA,textTitle);
fprintf(1, 'power Tc = %f', power_Tc);
%
%     textTitle = {['case:' num2str(nM) ' control:' num2str(nP) ' D\_adjust\_r, FGFR2 locus, 174 SNPs, control size 100, model group size 50'];'Null: M does not contain the individual'};
%     plot2hist(D_radjust_M0,D_radjust_MA,textTitle)
%
%     % textTitle = {['case' num2str(nM) 'control' num2str(nP) 'Statistic for SNP frequency, 2000 SNP, control size 100, model group size 50';'Null: M does not contain the individual']};
%     textTitle = {['case:' num2str(nM) 'control:' num2str(nP) ' D\_p ' fastafile];'Null: M does not contain the individual'};
%     plot2hist(D_p_M0,D_p_MA,textTitle);
%     textTitle = {'D\_rsign, FGFR2 locus, 174 SNPs, control size 100, model group size 50';'Null: M does not contain the individual'};
%     plot2hist(D_rsign_M0,D_rsign_MA,textTitle)
%     textTitle = {'D\_psign, FGFR2 locus, 174 SNPs, control 100, model group 50';'Null: M does not contain the individual'};
%     plot2hist(D_psign_M0,D_psign_MA,textTitle)
%     textTitle = {'D\_rcorrection, FGFR2 locus, 174 SNPs, control size 100, model group size 50';'Null: M does not contain the individual'};
%     plot2hist(D_rcorrection_M0,D_rcorrection_MA,textTitle)
%     textTitle = {'D\_pcorrection, FGFR2 locus, 174 SNPs, control 100, model group 50';'Null: M does not contain the individual'};
%     plot2hist(D_pcorrection_M0,D_pcorrection_MA,textTitle)
end



