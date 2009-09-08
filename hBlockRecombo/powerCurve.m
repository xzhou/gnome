function [] = powerCurve(allSeq, threshold, N)
% function powerCurve will test the power of Tr with Homer's attack.
% @caseSeq is the case hyplotype sequence
% @refSeq is reference hyplotype sequence
% @threshold filter r value that larger than threshold
% @resolution is the r value resolution
% @N is the N value
    
    nP = N;
    nM = N;
    
    Trials = 1000;  %the number of trilas 
    estimateSign = 0; % estimate the sign of r in the target sample from the reference
    removebySign = 0; % 1, remove the r values with the same sign in sample and reference.
                        % -1, remove the r values with different sign in sample
                        % and references
                        % else, do not remove r values by sign
    partialReverseSign = 0; % partially reverse the sign of the "removebySign" part; or totaly remove the "removebySign" part
    pseudocount = 0;    %remove MC bias
    precision = 10;
    
    if threshold > 0
        removeSmallr = 1;
    else
        removeSmallr = 0;
    end
    
    m = length(allSeq);
    n = length(allSeq(1).Sequence);
    
    seq4 = zeros(m, n);
    for i = 1:m
        seq4(i,:) = nt2int(allSeq(i).Sequence) - 1;
    end
    %allele mapping
    alleleMapping = getMajorAllele(seq4);
    
    %generate MC model from sequence
    try
        load('iMCmodel_SIM_1000x77.mat');
    catch ME
        iMCmodel = iMCmodelBuild(allSeq, pseudocount);
        save('iMCmodel_SIM_1000x77.mat', 'iMCmodel');
    end
    
    Len = size(iMCmodel.transition, 3) + 1;
    
    %sample one individual from model
    A = iMCgenerate(iMCmodel, 1);
    A = (A == alleleMapping) + 0;
    A1 = (2*A' - 1);
    A2 = (2*A'-1)*(2*A-1);
    
    % this is the test statistics each round
    D_r_M0 = zeros(Trials,1);
    D_r_MA = zeros(Trials,1);
    D_p_M0 = zeros(Trials,1);
    D_p_MA = zeros(Trials,1);
    
    %start multiple worker
    isOpen = matlabpool('size') > 0;
    if ~isOpen
        matlabpool 2;
    end
    
    countsmall = 0.0;
    
    %% repeat for 10000 times and calculate the distribution of r_M0 r_MA
    parfor i = 1:Trials
        %generate population from model
        sample_P = iMCgenerate(iMCmodel,nP);
        sample_M = iMCgenerate(iMCmodel,nM);

        sample_P = (sample_P == repmat(alleleMapping,nP,1)) + 0;
        sample_M = (sample_M == repmat(alleleMapping,nM,1)) + 0;

        all_r_P = corrcoef(sample_P);
        all_r_P(isnan(all_r_P))= 0; % deal with invariant sites

        all_r_M0 = corrcoef(sample_M);
        all_r_M0(isnan(all_r_M0))= 0; % deal with invariant sites

        all_p_P = sum(sample_P,1)/nP;
        all_p_M0 = sum(sample_M,1)/nM;
        
        %add the individual to the sample
        sample_M(end,:)=A;
        
        all_r_MA = corrcoef(sample_M);  %calculate r value
        all_r_MA(isnan(all_r_MA))= 0;   %deal with invariant sites
        all_p_MA = sum(sample_M,1)/nM;  %calculate single allele frequence

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
            sign_agree_M0P = sum(sum(sign(all_r_M0)==sign(all_r_P)))/numel(all_r_M0);
            sign_agree_MAP = sum(sum(sign(all_r_MA)==sign(all_r_P)))/numel(all_r_M0);
            sign_agree_MAM0 = sum(sum(sign(all_r_MA)==sign(all_r_M0)))/numel(all_r_M0);
            %replace sign here
            signMatrix = sign(all_r_P);
        end
        
        indx0 = [];
        indxA = []; 
        if removebySign == 1
            indx0 = sign(all_r_M0)==sign(all_r_P); %same as the reference
            indxA = sign(all_r_MA)==sign(all_r_P);  %same as the reference
        elseif removebySign == -1
            indx0 = sign(all_r_M0)==sign(all_r_P);  
            indxA = sign(all_r_MA)==sign(all_r_P);
        end
        
        if partialReverseSign == 1
            %partial ReverseSign
            all_r_M0(indx0) = -all_r_M0(indx0);
            all_r_MA(indxA) = -all_r_MA(indxA); 
        else
            %????
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

    
    %%normalize the Tr test
    D_r_M0 = D_r_M0/sqrt(Len*(Len-1)/2);
    D_r_MA = D_r_MA/sqrt(Len*(Len-1)/2);
    mean_r0 = mean(D_r_M0);
    mean_rA = mean(D_r_MA);
    var_r0 = var(D_r_M0);
    var_rA = var(D_r_MA);

    %%normalize the Tp test
    D_p_M0 = D_p_M0/sqrt(Len);
    D_p_MA = D_p_MA/sqrt(Len);

    countsmall = countsmall/Trials/((Len-1)*Len);

    % close all;
    % textTitle = {'Statistic for r, 2000 SNP, control size 100, model group size 50';'Null: M does not contain the individual'};
    
    %plot Tr
    
    textTitle = {['SIM 200x77 Power of Tr Threshold = ', num2str(threshold)]};
    power_Dr=plot2hist(D_r_M0,D_r_MA,textTitle);
    xlabel({['SIM 200x77 Power of Tr ' num2str(power_Dr)]});
    ylabel('Density of T_r');
    
    textTitle = {'SIM 200x77 Power of Tp Threshold = ', num2str(threshold)};
    power_Dp = plot2hist(D_p_M0,D_p_MA,textTitle);
    xlabel({['SIM 200x77 Power of Tp ' num2str(power_Dp)]});
    ylabel('Density of T_p');

end