function out = PowerCompare(fastafile, threshold, precision, N)
    %*P means reference
    % MA means sample
    % M0 means reference
    
    %cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77'
    cd 'D:\IUBResearch\Projects\Bioinfor\data\dist100x77';
    nP = 200;
    nM = 200;
    if nargin >=4
        nP = N;
        nM = N;
    end
    Trials = 1000;
    Len = 100;
    flag_usereal = 1;
    degFreedomFlag = 0; % estimate the degree of freedom of the haplotypes
    fastafile = '80SNP_CEU_sim_4000seq.fasta';
    % fastafile = 'chr10_FGFR2_200kb_phased_yri.fasta';
    frompairwise = 0; % compute the statistics from pairwise freq
    % fastafile = 'chr10_FGFR2_200kb_phased_jpt+chb.fasta';
    % fastafile =  'chr10_FGFR2_200kb_phased_CEU.fasta';
    % Tr_Tp = 1; % use Tr+Tp instead of Tr
    estimateSign = 1; % estimate the sign of r in the target sample from the reference
    removebySign = 0; % 1, remove the r values with the same sign in sample and reference.
                        % -1, remove the r values with different sign in sample
                        % and references
                        % else, do not remove r values by sign
    partialReverseSign = 0; % partially reverse the sign of the "removebySign" part; or totaly remove the "removebySign" part
    pseudocount = 0;
    
    if nargin < 2 || length(threshold)==0
        threshold = 0.1; %0.5, 0.2, 0.1, 0.01, 0.001
    end
    
    removeSmallr = 1;
    flag_twoPopulation = 0; % two models will be simulated for the references and the sample
    
    if nargin < 3 || length(precision)==0
        precision = 10; %1, 2, 3, 4, 5
    end

    %% get one model
    if ~ flag_usereal
    % Draw a model randomly
        nState = 2;
        iMCmodel = iMCmodelGen(Len, nState);
    else
    % or get from real data
        iMCmodel = iMCmodelBuild(fastafile, pseudocount);
        Len = size(iMCmodel.transition,3)+1;
        % nState = length(iMCmodel.initial);
    end


    % define allele as allele 1 for each of the SNP location, thus that SNP
    % with non 0/1 label can be converted to 0/1 label

    % allele mapping
    allele1 = iMCgenerate(iMCmodel,1);
    
    %???? why?
    if degFreedomFlag ==1
        d = zeros(200,1);
        for i = 1:200
            d(i) = sum(iMCgenerate(iMCmodel,1)==iMCgenerate(iMCmodel,1));
        end
        mean_d = mean(d);
        effictive_deg_freedom = (Len - mean_d)*2;
    end


    % Sample one individual form the model,
    A = iMCgenerate(iMCmodel,1);
    A = (A == allele1) + 0;
    A1 = (2*A'-1);
    A2 = (2*A'-1)*(2*A-1);  %sign

    % this is the test statistics each round
    D_r_M0 = zeros(Trials,1);
    D_r_MA = zeros(Trials,1);
    D_p_M0 = zeros(Trials,1);
    D_p_MA = zeros(Trials,1);
    
    
    %????
    countsmall = 0;
    if flag_twoPopulation ==1
        iMCmodel1 = iMCmodelBuild(iMCgenerate(iMCmodel,50),0);
        iMCmodel2 = iMCmodelBuild(iMCgenerate(iMCmodel,50),0);
    end
    i=1;
    
    %% repeat for 10000 times and calculate the distribution of r_M0 r_MA
    for i = 1: Trials
        if flag_twoPopulation ==1
            sample_P = iMCgenerate(iMCmodel1,nP);
            sample_M = iMCgenerate(iMCmodel2,nM);
        else
            sample_P = iMCgenerate(iMCmodel,nP);
            sample_M = iMCgenerate(iMCmodel,nM);
        end

        sample_P = (sample_P == repmat(allele1,nP,1)) + 0;
        sample_M = (sample_M == repmat(allele1,nM,1)) + 0;

        all_r_P = corrcoef(sample_P);
        all_r_P(isnan(all_r_P))= 0; % deal with invariant sites

        all_r_M0 = corrcoef(sample_M);
        all_r_M0(isnan(all_r_M0))= 0; % deal with invariant sites

        all_p_P = sum(sample_P,1)/nP;
        all_p_M0 = sum(sample_M,1)/nM;

        if frompairwise == 1
            [Tc_M0(i) Tc_log_M0(i) Tr_M0(i) TD_M0(i)] = Tc(sample_M, sample_P, A);
            sample_M(end,:)=A;
            [Tc_MA(i) Tc_log_MA(i) Tr_MA(i) TD_MA(i)] = Tc(sample_M, sample_P, A);
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
            all_r_M0(abs(all_r_M0)<threshold) = 0;
            all_r_MA(abs(all_r_MA)<threshold) = 0;
            all_r_P(abs(all_r_P)<threshold) = 0;
        end

        %calculate the pairwise test
        D_r_M0(i) = sum(sum((all_r_M0 - all_r_P).* A2))/2;
        D_r_MA(i) = sum(sum((all_r_MA - all_r_P).* A2))/2;

        % Test statistic for estimated frequency
        D_p_M0(i) = (all_p_M0 - all_p_P)* A1;
        D_p_MA(i) = (all_p_MA - all_p_P)* A1;
    end

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
    textTitle = {['case:' num2str(nM) ' control:' num2str(nP) ' Dr, FGFR2 locus, 174 SNPs, control size 100, model group size 50'];'Null: M does not contain the individual'};
    power_Dr=plot2hist(D_r_M0,D_r_MA,textTitle);
    xlabel('T_r');
    ylabel('Density of T_r');
    title('')

    textTitle = {['case:' num2str(nM) ' control:' num2str(nP) ' Dr, FGFR2 locus, 174 SNPs, control size 100, model group size 50'];'Null: M does not contain the individual'};
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

    textTitle = {['case:' num2str(nM) 'control:' num2str(nP) ' D\_p ' fastafile];'Null: M does not contain the individual'};
    power_Dp = plot2hist(D_p_M0,D_p_MA,textTitle);
    xlabel('T_p');
    ylabel('Density of T_p');
    title('')

    out = [threshold; countsmall; power_Dr; power_Dp; power_DrDp];
    textTitle = {['case:' num2str(nM) ' control:' num2str(nP) ' Tc: ' fastafile];'Null: M does not contain the individual'};
    plot2hist(Tc_M0,Tc_MA,textTitle)

    textTitle = {['case:' num2str(nM) ' control:' num2str(nP) ' Tc\_log: ' fastafile];'Null: M does not contain the individual'};
    plot2hist(Tc_log_M0,Tc_log_MA,textTitle)

    textTitle = {['case:' num2str(nM) ' control:' num2str(nP) ' D\_adjust\_r, FGFR2 locus, 174 SNPs, control size 100, model group size 50'];'Null: M does not contain the individual'};
    plot2hist(D_radjust_M0,D_radjust_MA,textTitle)

    % textTitle = {['case' num2str(nM) 'control' num2str(nP) 'Statistic for SNP frequency, 2000 SNP, control size 100, model group size 50';'Null: M does not contain the individual']};
    textTitle = {['case:' num2str(nM) 'control:' num2str(nP) ' D\_p ' fastafile];'Null: M does not contain the individual'};
    plot2hist(D_p_M0,D_p_MA,textTitle);
    textTitle = {'D\_rsign, FGFR2 locus, 174 SNPs, control size 100, model group size 50';'Null: M does not contain the individual'};
    plot2hist(D_rsign_M0,D_rsign_MA,textTitle)
    textTitle = {'D\_psign, FGFR2 locus, 174 SNPs, control 100, model group 50';'Null: M does not contain the individual'};
    plot2hist(D_psign_M0,D_psign_MA,textTitle)
    textTitle = {'D\_rcorrection, FGFR2 locus, 174 SNPs, control size 100, model group size 50';'Null: M does not contain the individual'};
    plot2hist(D_rcorrection_M0,D_rcorrection_MA,textTitle)
    textTitle = {'D\_pcorrection, FGFR2 locus, 174 SNPs, control 100, model group 50';'Null: M does not contain the individual'};
    plot2hist(D_pcorrection_M0,D_pcorrection_MA,textTitle)
end




