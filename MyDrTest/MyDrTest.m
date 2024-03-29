function [StatS, StatR, StatT, Truth] = MyDrTest(rev, replaceSign, maskMatrix, signMatrix)
% MYDRTEST is used to test the sign recovered by the distribution sign 
%   Method and compare the power of Tr test against 

    cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77'
    
    %% arguments processing
    precision = 4;
    Trials = 100;
    estimateStd = 0;
    Trialsonce = 1000;
    stdFromTest = 1; % esitimate variance from the Test values, with some power loss.
    muFromTest = 0; % estimate mu from Test values
    if estimateStd == 0
        Trials = Trialsonce;
    end

    if nargin == 0
        rev = 1;
        replaceSign = 0;
        maskMatrix = NaN;
        signMatrix = NaN;
    end

%% read haplotype sequence and get the sequence
    sample = 'SIM_100x77_smp.fasta';
    reference = 'SIM_100x77_ctl.fasta';
    test = 'SIM_100x77_smp.fasta';
    
    groupname = 'SIM';
    outgroupname = 'SIM';
    
    %reading sequence
    seqS = fastaread(sample);    
    seqR = fastaread(reference);
    seqT = fastaread(test);
    
    nSample = length(seqS);             % number of individuals from fasta file
    nS = nSample;
    Len = length(seqS(1).Sequence);      % number of snps in the fasta faile
    
    nR = length(seqR);               %the number of 
    %
    nTest = length(seqT);
    nT = length(seqT);
    
    int4S = zeros(nSample, Len);
    int4R = int4S;
    
    %convert sample and reference sequence to int
    for i = 1:length(seqS)
        int4S(i,:) = nt2int(seqS(i).Sequence) - 1;
        int4R(i,:) = nt2int(seqR(i).Sequence) - 1;
    end

    %convert 
    int4T = zeros(nT,Len); % test group
    for i = 1:length(seqT)
        int4T(i,:) = nt2int(seqT(i).Sequence) - 1;
    end
    

    %convert SNP nucleotide to 0/1
    % define allele as allele 1 for each of the SNP location, thus that SNP
    % with non 0/1 label can be converted to 0/1 label
    % int4* is the 4 integer labeling of SNPs ,and int2* is the 2 integer
    % labeling of SNPs
    allele1 = int4T(end,:);
    allele1 = getMajorAllele(int4R);
    int2S = (int4S == repmat(allele1,nS,1)) + 0; 
    int2R = (int4R == repmat(allele1,nR,1)) + 0;
    int2T = (int4T == repmat(allele1,nT,1)) + 0;
    
    %% calculate the r and p
    all_r_S = corrcoef(int2S);
    all_r_S(isnan(all_r_S)) = 0;%
    
    %replace the sign of the sample
    if replaceSign == 1
        %all_r_S = abs(all_r_S).*signMatrix;
    end
    

    all_r_R = corrcoef(int2R);
    all_r_R(isnan(all_r_S)) = 0;
    
    %initialize 
    StatS.Tr = zeros(nSample,1);
    StatS.std = StatS.Tr;
    StatS.p = StatS.Tr;
    
    StatR.Tr = zeros(nSample, 1);
    StatR.std = StatR.Tr;
    StatR.p = StatR.Tr;
    
    StatT.Tr = zeros(nSample, 1);
    StatT.std = StatT.Tr;
    StatT.p = StatT.Tr;
    
    %calculate the statistics for each individual
    if replaceSign == 1
        for i = 1:nSample
            StatS.Tr(i) = getTr(int2S(i,:), all_r_S, all_r_R, maskMatrix, signMatrix);
            StatR.Tr(i) = getTr(int2R(i,:), all_r_S, all_r_R, maskMatrix, signMatrix);
        end

        for i = 1:nTest
            StatT.Tr(i) = getTr(int2T(i,:), all_r_S, all_r_R, maskMatrix, signMatrix);
        end
    else
        for i = 1:nSample
            StatS.Tr(i) = getTr(int2S(i,:), all_r_S, all_r_R);
            StatR.Tr(i) = getTr(int2R(i,:), all_r_S, all_r_R);
        end

        for i = 1:nTest
            StatT.Tr(i) = getTr(int2T(i,:), all_r_S, all_r_R);
        end
    end
    
    %????
    StatS.Tr = StatS.Tr/sqrt(Len*(Len-1)/2);
    StatR.Tr = StatR.Tr/sqrt(Len*(Len-1)/2);
    StatT.Tr = StatT.Tr/sqrt(Len*(Len-1)/2);
    %calculate p
    
    %% read snp.plotter data and replace the sign

    %%plot
    Trall = [StatS.Tr; StatR.Tr; StatT.Tr];
    StdAll = [StatS.std; StatR.std; StatT.std];
    pall = [StatS.p; StatR.p; StatT.std];
    
    indexS = 1:nS;
    indexR = (nS+1):(nS+nR);
    indexT = (nS+nR+1):(nS+nR+nT);
    
    figure;
    hold on;
    
    nspace = 8;
    plot(indexS, Trall(indexS), '.r');
    plot(indexR, Trall(indexR), '.g');
    %plot(indexT, Trall(indexT), '.k');
    a = axis;
    xlabel('index of individuals');
    ylabel('T_r Values');
    legend({['sample'] ['Ref']});
    saveas(gcf, ['Tr_', groupname '_vs_' outgroupname '_' '.fig']);
    saveas(gcf, ['Tr_', groupname '_vs_' outgroupname '_' '.png']);
    
end


function [std_Dr, mean_Dr] = getstdTr(int4, iMCmodel, refhaplotype, nS, nR, Trials, precision)
% UNCLEAR, I guess it is get the standard deviation from reference
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

        %sample_S(end,:)=A;
        all_r_R = corrcoef(sample_R);
        all_r_R(isnan(all_r_R))= 0; % deal with invariant sites
        all_r_S = corrcoef(sample_S);
        all_r_S(isnan(all_r_S))= 0; % deal with invariant sites

        % Test statistic for r
        D_r_S(i) = sum(sum((all_r_S - all_r_R).* A2))/2;
    end
    D_r_S = D_r_S/sqrt(Len*(Len-1)/2);
    % figure
    % hist(D_r_S);
    mean_Dr = mean(D_r_S);
    std_Dr = sqrt(var(D_r_S));
    return
end

function Tr = getTr(Y, r_S, r_R, maskMatrix, signMatrix)
    if nargin == 3
        A2 = (2*Y'-1)*(2*Y-1);
        Tr = sum(sum((r_S - r_R).* A2))/2;
    else
        %r_S = abs(r_S).*signMatrix; %replace the sign
        r_S = r_S.*maskMatrix;
        r_R = r_R.*maskMatrix;
        A2 = (2*Y'-1)*(2*Y-1);
        Tr = sum(sum((r_S - r_R).* A2))/2;
    end
end