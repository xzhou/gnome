function [StatS, StatR, StatT, graph] = powerTest(caseSeq, refSeq, alleleMapping, caseR, maskMatrix, signMatrix)
% MYDRTEST is used to test the sign recovered by the distribution sign 
%   Method and compare the power of Tr test against 
%   @caseSeq    [0 1 2 3] to represent case [A C T G];
%   @refSeq     reference sequence
%   @alleleMapping  the major allele mapping
%   @caseR          caseR can come from R dump data
%   @maskMatrix     mask the r position that we don't want to count 0 means
%                   ignore
%   @signMatrix     use the recoverd sign to replace the sing in case R

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
    
    if nargin == 2
        alleleMapping = getMajorAllele(refSeq);
        refR = calcR(refSeq, alleleMapping);
        caseR = calcR(caseSeq, alleleMapping);
        replaceSign = 0;
        maskMatrix = zeros(size(caseR));
        signMatrix = zeros(size(caseR));
    elseif nargin == 3
        refR = calcR(refSeq, alleleMapping);
        caseR = calcR(caseSeq, alleleMapping);
        replaceSign = 0;
        maskMatrix = zeros(size(caseR));
        signMatrix = zeros(size(caseR));
    elseif nargin == 4
        refR = calcR(refSeq, alleleMapping);
        replaceSign = 0;
        maskMatrix = zeros(size(caseR));
        signMatrix = zeros(size(caseR));
    elseif nargin == 6
        replaceSign = 1;
        refR = calcR(refSeq, alleleMapping);
    end
    
    
    %%statistics
    [m n] = size(caseR);
    nLarger01 = sum(sum(caseR>=0.1));
    
    
    
    %%convert the sequence to 0/1 representation
    nS = length(caseSeq);
    nR = length(refSeq);
    Len = length(caseSeq(1,:));
    int2S = (caseSeq == repmat(alleleMapping, nS, 1)) + 0;
    int2R = (refSeq == repmat(alleleMapping, nR, 1)) + 0;
    
    %%initialize 
    StatS.Tr = zeros(nS,1);
    StatS.std = StatS.Tr;
    StatS.p = StatS.Tr;
    
    StatR.Tr = zeros(nR, 1);
    StatR.std = StatR.Tr;
    StatR.p = StatR.Tr;
    
    %calculate the statistics for each individual
    detailCase = zeros(Len);
    detailRef = zeros(Len);
    
    if replaceSign == 1
        for i = 1:nS
            [StatS.Tr(i) detail] = getTr(int2S(i,:), caseR, refR, maskMatrix, signMatrix);
%             detail(detail>0) = 1;
%             detail(detail<=0) = 0;
            detailCase = detailCase + detail;
            [StatR.Tr(i) detail] = getTr(int2R(i,:), caseR, refR, maskMatrix, signMatrix);
%             detail(detail>0) = 1;
%             detail(detail<=0) = 0;
            detailRef = detailRef + detail;
        end
    else
        for i = 1:nS
            [StatS.Tr(i) detail] = getTr(int2S(i,:), caseR, refR);
            detailCase = detailCase + detail;
            [StatR.Tr(i) detail] = getTr(int2R(i,:), caseR, refR);
            detailRef = detailRef + detail;
        end
    end
    
    signDiff = sign(caseR).*sign(refR) + 0;
    signDiff(~maskMatrix) = 0;
    
%     plotPowerDist(detailCase, signDiff);
%     saveas(gcf, 'casePowerDist.pdf');
%     plotPowerDist(detailCase);
%     saveas(gcf, 'casePowerDistNoSign.pdf');
%     plotPowerDist(detailRef, signDiff);
%     saveas(gcf, 'refPowerDist.pdf');
%     plotPowerDist(detailRef);
%     saveas(gcf, 'refPowerDistNoSign.pdf');
    
    plotPowerDist(detailCase - detailRef);
    saveas(gcf, 'powerDistNoSign 0.1.pdf');
    plotPowerDist(detailCase - detailRef, signDiff);
    saveas(gcf, 'powerDistWithSign 0.1.pdf');
    
    %????
    StatS.Tr = StatS.Tr/sqrt(Len*(Len-1)/2);
    StatR.Tr = StatR.Tr/sqrt(Len*(Len-1)/2);

    %%plot
    Trall = [StatS.Tr; StatR.Tr];
    
    indexS = 1:nS;
    indexR = (nS+1):(nS+nR);
    
    figure;
    hold on;
    
    nspace = 8;
    plot(indexS, Trall(indexS), '.r');
    plot(indexR, Trall(indexR), '.g');
    a = axis;
    xlabel('index of individuals');
    ylabel('T_r Values');
    legend({['sample'] ['Ref']});
    %saveas(gcf, ['Tr_', groupname '_vs_' outgroupname '_' '.fig']);
    %saveas(gcf, ['Tr_', groupname '_vs_' outgroupname '_' '.png']);
    graph = gcf;
end


%calculate the standard deviation for each Tr value, long time operation
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

function [Tr detail] = getTr(Y, r_S, r_R, maskMatrix, signMatrix)
    if nargin == 3
        A2 = (2*Y'-1)*(2*Y-1);
        Tr = sum(sum((r_S - r_R).* A2))/2;
        detail = (r_S - r_R).*A2;
    else
        %r_S = abs(r_S).*signMatrix; %replace the sign
        r_S = r_S.*maskMatrix;
        r_R = r_R.*maskMatrix;
        A2 = (2*Y'-1)*(2*Y-1);
        detail = (r_S - r_R).*A2;
        Tr = sum(sum((r_S - r_R).* A2))/2;
    end
end