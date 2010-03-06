%% do relaxed hyplotype block recombination for real data. 
%We will try to relax the consistency model for free recombination between
%any blocks to see if we can recover more signs;

%change_env();
%start computing slaves

startParallel();

blocks = [1, 24; 25, 65; 66, 81]; 

try
    cd 'D:\IUBResearch\Projects\Bioinfor\data\1500';
    disp 'WINDOWS'
catch e
    try 
        %cd '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/TPED';
        cd '~/GProject/1500/';
        disp 'LINUX'
    catch e
    end
end

[nBlock tmp] = size(blocks);

delete('hbrecombo.log');
diary hbrecombo.log;

rawFastaData = fastaread('newAffx.fasta');
%rawFastaData = fastaread('Affx_gt_58C_Chiamo_07.tped.fasta');

bigRepeat = 20;
averagePreSign = 0;
averagePostSign = 0;
averageEstiSign = 0;
averagePostEstiSign = 0;

averagePre95 = 0;
averagePost95 = 0;
averageCorrect95 = 0;
averageEsti95 = 0;
averageHomer95 = 0;

for iBigRepeat = 1: bigRepeat
    fprintf(1, '\n*** trial %d ***\n', iBigRepeat);
    
    
%      for 77 SNPs randomly select case and reference
    [caseSeq4 refSeq4] = randomSelect(rawFastaData, 250, 0);

%     caseFasta = fastaread('chr10_FGFR2_200kb_phased_CEU.fasta');
%     refFasta = fastaread('chr10_FGFR2_200kb_phased_yri.fasta');
%     
    nS = length(caseSeq4(:,1));
    Len = length(caseSeq4(1,:));
%     caseSeq4 = zeros(nS, Len);
%     refSeq4 = zeros(nS,Len);
%     
%     for i = 1 :nS
%         caseSeq4(i, :) = nt2int(caseFasta(i).Sequence) -1;
%         refSeq4(i,:) = nt2int(refFasta(i).Sequence) -1;
%     end


    %% test the block structure
    caseBlockFreqInfo = cell(nBlock, 1);

    parfor i = 1:nBlock
        caseBlockFreqInfo{i,1} = getBlockFreq(caseSeq4, blocks(i,:));
    end

    refBlockFreqInfo = cell(nBlock, 1);
    parfor i = 1:nBlock
        refBlockFreqInfo{i,1} = getBlockFreq(refSeq4, blocks(i,:));
    end
    

    matchedCase = blockCheck(caseBlockFreqInfo, refBlockFreqInfo, blocks);
    refMatchedCase = blockCheck(refBlockFreqInfo, caseBlockFreqInfo, blocks);
    
    

    for i = 1:nBlock
        [tempRows tempCols] = size(caseBlockFreqInfo{i,1});
        fprintf (1, '*********Block%d********* \n ', i);
        fprintf (1, 'Case    Hyplotype Blocks:%d \n ', tempRows);
        [tempRows tempCols] = size(refBlockFreqInfo{i,1});
        fprintf (1, 'Ref     Hyplotype Blocks:%d \n ', tempRows);
        fprintf (1, 'Shared  Hyplotype Blocks:%d \n ', sum(matchedCase{i,1}(end,:)~=0));
    end
    
    save('refBlockFreq.mat', 'refBlockFreqInfo');
    save('caseBlockFreq.mat', 'caseBlockFreqInfo');
    save('caseMatch.mat', 'matchedCase');
    save('refMatchedCase.mat', 'refMatchedCase');

    % %% Construct current sequence accorting to CASE frequence
    % 
    % afterInnerRecomboSeq = zeros(nS, Len);
    % tempRefMatchedCase = refMatchedCase;
    % n = 1;
    % for i = 1 : length(refMatchedCase)
    %     tempRefMatchedCase{i, 1}(end-1:end, :) = [];
    %     for j = 1: length(tempRefMatchedCase{i,1})
    %         if tempRefMatchedCase{i,1}(j, end) ~= 0
    %             for k = 1 : tempRefMatchedCase{i,1}(j, end)
    %             afterInnerRecomboSeq(n, blocks(i,1):blocks(i,2)) = tempRefMatchedCase{i,1}(j, 1:end-2);
    %             n = n+1;
    %             end
    %         end
    %     end
    %     n= 1;
    % end


    %% plot initialization

    %index for plotting
    index1 = [1: nS];
    index2 = [nS+1: nS*2];

    alleleMapping = getMajorAllele(refSeq4);

    int2S = (caseSeq4 == repmat(alleleMapping,nS,1)) + 0;
    int2R = (refSeq4 == repmat(alleleMapping,nS,1)) + 0;

     %% Test the power of Homer's test
    singleFreS = sum(int2S, 1)/nS;
    singleFreR = sum(int2R, 1)/nS;
    StatS.Tp = zeros(nS, 1);
    StatR.Tp = zeros(nS, 1);

    for i = 1: nS
        StatS.Tp(i) = (singleFreS - singleFreR)*(2*int2S(i,:)' - 1);
        StatR.Tp(i) = (singleFreS - singleFreR)*(2*int2R(i,:)' - 1);
    end
    sortHomerStatR = sort(StatR.Tp);
    homerAbove95S = sum(StatS.Tp>sortHomerStatR(int16(nS*0.95)));
    h=figure;
    hold on;
    plot(index1, StatS.Tp, '.r');
    plot(index2, StatR.Tp, '.g');
    legend({'case' 'ref'});
    plot(ones(2*nS).*sortHomerStatR(int16(nS*0.95)));
    xlabel('individual index');
    ylabel('T_r value');
    title('Homer Test');
    filename = strcat('Trial', num2str(iBigRepeat),'HomerTest.pdf');
    print(h,'-dpdf',filename);

    caseGenoSeq = haplotype2genotype(caseSeq4, alleleMapping);
    [targetR pA counts] = estimateR(caseGenoSeq);
	targetR = calcR(caseSeq4, alleleMapping);
	
    nanAdjust = sum(sum(isnan(targetR)));
	%targetR(isnan(targetR))=0;
    realSingleAlleleFreq = GnomeCalculator.getSingleAlleleFreq(caseSeq4, alleleMapping);
    %round to precision 1e-4
    targetR = fix(targetR.*10000)./10000;
    realTargetR = calcR(caseSeq4, alleleMapping);
    
    refR = calcR(refSeq4, alleleMapping);
    
	[refR2 refC00 refC01 refC10 refC11] = calcPairwiseFreq(refSeq4, alleleMapping);

    preTargetR = abs(targetR).*sign(refR);

    preStatS.Tr = zeros(nS, 1);
    preStatR.Tr = zeros(nS, 1);

    postStatS.Tr = zeros(nS, 1);
    postStatR.Tr = zeros(nS, 1);

    correctStatS.Tr = zeros(nS, 1);
    correctStatR.Tr = zeros(nS, 1);
    
    estimateValueStatS.Tr = zeros(nS, 1);
    estimateValueStatR.Tr = zeros(nS, 1);
    
    
    %For caculating the Tr with the correct Sign and Value from EstimateR
    estimateValueR = abs(targetR).*sign(realTargetR);
    for i = 1:nS
        estimateValueStatS.Tr(i) = getTr(int2S(i,:), estimateValueR, refR);
        estimateValueStatR.Tr(i) = getTr(int2R(i,:), estimateValueR, refR);
    end
    estimateValueStatS.Tr = estimateValueStatS.Tr/sqrt(Len*(Len-1)/2);
    estimateValueStatR.Tr = estimateValueStatR.Tr/sqrt(Len*(Len-1)/2);


    sortestimateValueStatR = sort(estimateValueStatR.Tr);
    estimateValueAbove95S = sum(estimateValueStatS.Tr>sortestimateValueStatR(int16(nS*0.95)));


    %For caculating the Tr with correct Sign and correct R value
    for i = 1:nS
        correctStatS.Tr(i) = getTr(int2S(i,:), realTargetR, refR);
        correctStatR.Tr(i) = getTr(int2R(i,:), realTargetR, refR);
    end
    correctStatS.Tr = correctStatS.Tr/sqrt(Len*(Len-1)/2);
    correctStatR.Tr = correctStatR.Tr/sqrt(Len*(Len-1)/2);


    sortCorrectStatR = sort(correctStatR.Tr);
    correctAbove95S = sum(correctStatS.Tr>sortCorrectStatR(int16(nS*0.95)));

    %For caculating the Tr befor recombination (with Sign from Ref and
    %estimate value
    for i = 1:nS
        preStatS.Tr(i) = getTr(int2S(i,:), preTargetR, refR);
        preStatR.Tr(i) = getTr(int2R(i,:), preTargetR, refR);
    end
    preStatS.Tr = preStatS.Tr/sqrt(Len*(Len-1)/2); 
    preStatR.Tr = preStatR.Tr/sqrt(Len*(Len-1)/2);

    sortPreStatR = sort(preStatR.Tr);
    preAbove95S = sum(preStatS.Tr>sortPreStatR(int16(nS*0.95)));

    %finalTargetR = preTargetR;


    %Some problems happens here about blocks??

    %for blocks larger than 3
    blocks = [blocks blocks(:,2) - blocks(:,1) + 1];
    [m n] = size(blocks);
    %add learning marker
    blocks = [blocks, zeros(m, 1)];
    %blocks = blocks(blocks(:,3)>=3, :);


    [m n] = size(blocks);
    finalResult = zeros(m, m);

    trials = 16;

    f = fopen('result.txt', 'w');

    %Save all the sign matrix
    bufferMatrix = zeros(Len, Len, trials);

    %% inner block learning
    %currentSeq = innerBlockDriver(caseSeq4, refSeq4, blocks);
    currentSeq = refSeq4;
    
    currentBlockFreqInfo = cell(nBlock, 1);

    parfor i = 1:nBlock
        currentBlockFreqInfo{i,1} = getBlockFreq(currentSeq, blocks(i,:));
    end
    
    caseMatchedCurrent = blockCheck(caseBlockFreqInfo, currentBlockFreqInfo, blocks);
    currentMatchedCase = blockCheck(currentBlockFreqInfo, caseBlockFreqInfo, blocks);
    

    newCurrentSeq = adjustStartPoint(currentSeq, refSeq4, blocks);

    finalTargetR = calcR(newCurrentSeq, alleleMapping);
    
	finalTargetR = abs(targetR).*sign(finalTargetR);
    
    %signMatrix = zeros(size(caseSeq4));
    %signMatrix = sign(calcR(currentSeq));
    for i = 1:(m-1)
        for j = i+1:m
            blockRate = zeros(trials, 2);      
            block1 = blocks(i,:);
            block2 = blocks(j,:);
            if(block1(1,3) >= block2(1,3))
                block1(1,4) = 1;%??
                %currentSeq = shuffleNewBlock(currentSeq, block2);
                for t = 1:trials
                    [finalSeq finalR finalSignRate finalQual blockMask] = newHBRecombo(targetR, caseSeq4, refSeq4, newCurrentSeq, block1, block2, alleleMapping, 0);
                    bufferMatrix(:,:,t) = sign(finalR.*blockMask);
                    blockRate(t,:) = [finalSignRate finalQual];
                    if finalQual == 0 && finalSignRate ~= 1.0
                        pause
                    end
                end
            else
                block2(1,4) = 1;%??
                %currentSeq = shuffleNewBlock(currentSeq, block1);
                for t = 1:trials
                    [finalSeq finalR finalSignRate finalQual blockMask] = newHBRecombo(targetR, caseSeq4, refSeq4, newCurrentSeq, block2, block1, alleleMapping, 0);
                    bufferMatrix(:,:,t) = sign(finalR.*blockMask);
                    blockRate(t,:) = [finalSignRate finalQual];
                    if finalQual == 0 && finalSignRate ~= 1.0
                        pause
                    end
                end
            end


            %for max recover rate
            [maxVal, maxIdx] = max(blockRate(:,1));

            %optimal
            [minQ, minIdx] = min(blockRate(:,2));
            maxQ = min(blockRate(:,2));

            %Apply signs to these cross blocks
            finalTargetR = finalTargetR.*(bufferMatrix(:,:,minIdx)==0) + abs(finalTargetR).*(bufferMatrix(:,:,minIdx)~=0).*bufferMatrix(:,:,minIdx);

            finalResult(i,j) = maxVal;

            fprintf(1, 'optimal r dist = %f, signRate = %f\n', blockRate(minIdx,2),blockRate(minIdx,1));
    %         if maxQ == blockRate(maxIdx, 2)
    %             fprintf(f, '[%d-%d]x[%d-%d]\t = %f\t YES\n', block1(1,1), block1(1,2), block2(1,1), block2(1,2), finalResult(i,j));
    %             fprintf(1, '[%d-%d]x[%d-%d]\t = %f\t YES\n', block1(1,1), block1(1,2), block2(1,1), block2(1,2), finalResult(i,j));
    %         else
    %             fprintf(f, '[%d-%d]x[%d-%d]\t = %f\t NO\n', block1(1,1), block1(1,2), block2(1,1), block2(1,2), finalResult(i,j));
    %             fprintf(1, '[%d-%d]x[%d-%d]\t = %f\t NO\n', block1(1,1), block1(1,2), block2(1,1), block2(1,2), finalResult(i,j));
    %         end
        end
    end

    %% For caculating the Tr after recombination
    for i = 1:nS
        postStatS.Tr(i) = getTr(int2S(i,:), finalTargetR, refR);
        postStatR.Tr(i) = getTr(int2R(i,:), finalTargetR, refR);
    end

    postStatS.Tr = postStatS.Tr/sqrt(Len*(Len-1)/2);
    postStatR.Tr = postStatR.Tr/sqrt(Len*(Len-1)/2);

    sortPostStatR = sort(postStatR.Tr);
    postAbove95S = sum(postStatS.Tr>sortPostStatR(int16(nS*0.95)));

    plotScatter(caseSeq4, refSeq4, preTargetR, refR, 'WithSignsfromREF', iBigRepeat);
    plotScatter(caseSeq4, refSeq4, finalTargetR, refR, 'AfterSignRecovery', iBigRepeat);
    plotScatter(caseSeq4, refSeq4, targetR, refR, 'WithCorrectSigns', iBigRepeat);
	
	[preSignRate test1 test2] = SignRate(realTargetR, refR);
	[postSignRate test1 test2] = SignRate(realTargetR, finalTargetR);
	[estiSignRate test1 test2] = SignRate(realTargetR, targetR);
	[postToEstiSignRate test1 test2] = SignRate(finalTargetR, targetR);
    
    % preSignRate = (sum(sum(sign(realTargetR)==sign(refR)))-nanAdjust)/(Len*(Len-1)-nanAdjust);
    % postSignRate = (sum(sum(sign(realTargetR)==sign(finalTargetR)))-nanAdjust)/(Len*(Len-1)-nanAdjust);
	% estiSignRate = (sum(sum(sign(realTargetR)==sign(targetR)))-nanAdjust)/(Len*(Len-1)-nanAdjust);
	% postToEstiSignRate = (sum(sum(sign(finalTargetR)==sign(targetR)))-nanAdjust)/(Len*(Len-1)-nanAdjust);
    
    %preSignRate = sum(sum(sign(targetRealR)==sign(refR)))/Len/(Len-1);
    %postSignRate = sum(sum(sign(targetRealR)==sign(finalTargetR)))/Len/(Len-1);
	
	averagePreSign = averagePreSign + preSignRate;
	averagePostSign = averagePostSign + postSignRate;
	averageEstiSign = averageEstiSign + estiSignRate;
	averagePostEstiSign = averagePostEstiSign + postToEstiSignRate;
	
	averagePre95 = averagePre95 + preAbove95S;
	averagePost95 = averagePost95 + postAbove95S;
	averageCorrect95 = averageCorrect95 + correctAbove95S;
	averageEsti95 = averageEsti95 + estimateValueAbove95S;
	averageHomer95 = averageHomer95 + homerAbove95S;
    
    fprintf (1, 'PreSignRate = %f\n ', preSignRate);
    fprintf (1, 'PostSignRate = %f\n ', postSignRate);
	fprintf (1, 'EstiSignRate = %f\n ', estiSignRate);
	fprintf (1, 'PostToEstiSignRate = %f\n ', postToEstiSignRate);
	
    fprintf (1, '**************************************\n');
    fprintf (1, 'PreAbove95s             %d\n ', preAbove95S);
    fprintf (1, 'PostAbove95s            %d\n ', postAbove95S);
    fprintf (1, 'CorrectAbove95s         %d\n ', correctAbove95S);
    fprintf (1, 'CorrectSignEstimateRAbove95s            %d\n ', estimateValueAbove95S);
    fprintf (1, 'HomerAbove95s           %d\n ', homerAbove95S);
%     fprintf (1, '**************************************\n');
%     fprintf (1, 'PreCaseTr>0         %d\n ', sum(preStatS.Tr>0));
%     fprintf (1, 'PostCaseTr>0        %d\n ', sum(postStatS.Tr>0));
%     fprintf (1, 'PreReferenceTr<0    %d\n ', sum(preStatR.Tr<0));
%     fprintf (1, 'PostReferenceTr<0   %d\n ', sum(postStatR.Tr<0));
%     fprintf (1, 'CorrectCaseTr>0     %d\n ', sum(correctStatS.Tr>0));
%     fprintf (1, 'CorrectReference<0  %d\n ', sum(correctStatR.Tr<0));
end


    fprintf (1, '**************************************\n');
	fprintf (1, 'Average Test results for %d bigRepeat and %d trials\n' , bigRepeat, trials);
	fprintf (1, '**************************************\n');
	fprintf (1, 'AveragePreSignRate = %f\n ', averagePreSign/bigRepeat);
    fprintf (1, 'AveragePostSignRate = %f\n ', averagePostSign/bigRepeat);
	fprintf (1, 'AverageEstiSignRate = %f\n ', averageEstiSign/bigRepeat);
	fprintf (1, 'AveragePostToEstiSignRate = %f\n ', averagePostEstiSign/bigRepeat);
	fprintf (1, '**************************************\n');
    fprintf (1, 'AveragePreAbove95s             %f\n ', averagePre95/bigRepeat);
    fprintf (1, 'AveragePostAbove95s            %f\n ', averagePost95/bigRepeat);
    fprintf (1, 'AverageCorrectAbove95s         %f\n ', averageCorrect95/bigRepeat);
    fprintf (1, 'AverageCorrectSignEstimateRAbove95s            %f\n ', averageEsti95/bigRepeat);
    fprintf (1, 'AverageHomerAbove95s           %f\n ', averageHomer95/bigRepeat);

fclose(f);
diary off;

save;
save('finalResult.mat', 'finalResult');
save('caseSeq4.mat', 'caseSeq4');
save('refSeq4.mat', 'refSeq4');
save('finalTargetR.mat', 'finalTargetR');
