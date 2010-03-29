%% do relaxed hyplotype block recombination for real data.
%We will try to relax the consistency model for free recombination between
%any blocks to see if we can recover more signs;

%change_env();
%start computing slaves

startParallel();

%blocks = [1, 24; 25, 65; 66, 81];
%blocks = [1, 24];

%We may do some improvement here by dividing the seq into more blocks
%blocks = [1, 24; 25, 65; 66, 81; 82, 121; 122, 136; 137, 156; 157, 176; 177, 204]; 
%7 Blocks
%blocks = [1, 24; 25, 56; 57, 81; 82, 121; 122, 136; 137, 156; 157, 204]; 
%5 Blocks
blocks = [1, 24; 25, 81; 82, 121; 122, 136; 137, 204]; 
%5Blocks with 176 snps
% blocks = [1, 24; 25, 81; 82, 121; 122, 136; 137, 176];
%4 Blocks with 176 snps
% blocks = [1, 24; 25, 81; 82, 121; 122, 176]; 
 

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

%rawFastaData = fastaread('new200snp.fasta');
rawFastaData = fastaread('600Snp.fasta');
%rawFastaData = fastaread('Affx_gt_58C_Chiamo_07.tped.fasta');

bigRepeat = 50;

trials = 12;

averagePreSign = 0;
averagePreEstiSign = 0;
averagePostSign = 0;
averageEstiSign = 0;
averagePostEstiSign = 0;

averagePre95 = 0;
averagePost95 = 0;
averageCorrect95 = 0;
averageEsti95 = 0;
averageHomer95 = 0;
averageRealEsti95 = 0;

averagePre99 = 0;
averagePost99 = 0;
averageCorrect99 = 0;
averageEsti99 = 0;
averageHomer99 = 0;
averageRealEsti99 = 0;


for iBigRepeat = 1: bigRepeat
    fprintf(1, '\n*** trial %d ***\n', iBigRepeat);
    
    %% Get the samples from fasta data
    [nS ~] = size(rawFastaData);
    for i = 1:nS
    int4All(i,:) = nt2int(rawFastaData(i).Sequence) - 1;
    end
    
    alleleMapping = getMajorAllele(int4All);
    allSingleAlleleFreq = GnomeCalculator.getSingleAlleleFreq(int4All, alleleMapping);

    int4Unique = unique(int4All, 'rows');
    
    %[caseSeq4 refSeq4 testSeq4] = randomSelectTest(int4Unique, 100);
    
    [caseSeq4 refSeq4 testSeq4] = randomSelectTestNew(rawFastaData, 200);
    
%     caseSeq4(:,177:end)=[];
%     refSeq4(:,177:end)=[];
%     testSeq4(:,177:end)=[];
    
    nS = length(caseSeq4(:,1));
    Len = length(caseSeq4(1,:));
    
    
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
        fprintf (1, 'Shared  Hyplotype Blocks:%d \n ', sum(matchedCase{i,1}(:,end)~=0));
    end
    
    save('refBlockFreq.mat', 'refBlockFreqInfo');
    save('caseBlockFreq.mat', 'caseBlockFreqInfo');
    save('caseMatch.mat', 'matchedCase');
    save('refMatchedCase.mat', 'refMatchedCase');
    
    
    %% plot initialization
    
    %index for plotting
    index1 = [1: nS/2];
    index2 = [nS/2+1: nS];
    index3 = [nS+1:nS/2*3];
    
    alleleMapping = getMajorAllele(refSeq4);
    
    int2S = (caseSeq4 == repmat(alleleMapping,nS,1)) + 0;
    int2R = (refSeq4 == repmat(alleleMapping,nS,1)) + 0;
    int2T = (testSeq4 == repmat(alleleMapping,nS,1)) + 0;
    
    %% Test the power of Homer's test
    %Check the correctness of Homer's test
    singleFreS = sum(int2S, 1)/nS;
    singleFreR = sum(int2R, 1)/nS;
    singleFreT = sum(int2T, 1)/nS;
    
    alleleFreS = zeros(nS/2, Len);
    alleleFreR = zeros(nS/2, Len);
    alleleFreT = zeros(nS/2, Len);
    
    StatS.Tp = zeros(nS/2, 1);
    StatR.Tp = zeros(nS/2, 1);
    StatT.Tp = zeros(nS/2, 1);
    
    
    
    for i = 1: nS/2
        alleleFreS(i,:) = (int2S(i,:)+int2S(i+1,:))/2;
        alleleFreR(i,:) = (int2R(i,:)+int2R(i+1,:))/2;
        alleleFreT(i,:) = (int2T(i,:)+int2T(i+1,:))/2;
        
        StatS.Tp(i) = sum(abs(alleleFreS(i,:)-singleFreR)-abs(alleleFreS(i,:)-singleFreS));
        StatR.Tp(i) = sum(abs(alleleFreR(i,:)-singleFreR)-abs(alleleFreR(i,:)-singleFreS));
        StatT.Tp(i) = sum(abs(alleleFreT(i,:)-singleFreR)-abs(alleleFreT(i,:)-singleFreS));
%         StatS.Tp(i) = (singleFreS - singleFreR)*(2*int2S(i,:)' - 1);
%         StatR.Tp(i) = (singleFreS - singleFreR)*(2*int2R(i,:)' - 1);
%         StatT.Tp(i) = (singleFreS - singleFreR)*(2*int2T(i,:)' - 1);
    end
    
    %fprintf(1, 'Homer average Case:%f; Ref:%f; Test:%f\n', mean(StatS.Tp), mean(StatR.Tp), mean(StatT.Tp));
    
    sortHomerStatT = sort(StatT.Tp);
    homerAbove95S = sum(StatS.Tp>sortHomerStatT(int16(nS/2*0.95)));
    homerAbove99S = sum(StatS.Tp>sortHomerStatT(int16(nS/2*0.99)));
    h=figure;
    hold on;
    plot(index1, StatS.Tp, '.r');
    plot(index2, StatR.Tp, '.g');
    plot(index3, StatT.Tp, '.b');
    legend({'case' 'ref' 'test'});
    plot(ones(3*nS/2).*sortHomerStatT(int16(nS/2*0.95)));
    plot(ones(3*nS/2).*sortHomerStatT(int16(nS/2*0.99)));
    xlabel('individual index');
    ylabel('T_r value');
    title('Homer Test');
    filename = strcat('Trial', num2str(iBigRepeat),'HomerTest.pdf');
    print(h,'-dpdf',filename);

    caseGenoSeq = haplotype2genotype(caseSeq4, alleleMapping);
    %[targetR pA counts] = estimateR(caseGenoSeq);
	targetR = calcR(caseSeq4, alleleMapping);
	
    nanAdjust = sum(sum(isnan(targetR)));
	targetR(isnan(targetR))=0;
    realSingleAlleleFreq = GnomeCalculator.getSingleAlleleFreq(caseSeq4, alleleMapping);
    %round to precision 1e-4
    targetR = fix(targetR.*10000)./10000;
    realTargetR = calcR(caseSeq4, alleleMapping);
    
    refR = calcR(refSeq4, alleleMapping);
    
	[refR2 refC00 refC01 refC10 refC11] = calcPairwiseFreq(refSeq4, alleleMapping);

    preTargetR = abs(targetR).*sign(refR);
    
    preStatS.Tr = zeros(nS, 1);
    preStatR.Tr = zeros(nS, 1);
    preStatT.Tr = zeros(nS, 1);
    
    postStatS.Tr = zeros(nS, 1);
    postStatR.Tr = zeros(nS, 1);
    postStatT.Tr = zeros(nS, 1);
    
    correctStatS.Tr = zeros(nS, 1);
    correctStatR.Tr = zeros(nS, 1);
    correctStatT.Tr = zeros(nS, 1);  
    
    estimateValueStatS.Tr = zeros(nS, 1);
    estimateValueStatR.Tr = zeros(nS, 1);
    estimateValueStatT.Tr = zeros(nS, 1);
    
    estimateStatS.Tr = zeros(nS, 1);
    estimateStatR.Tr = zeros(nS, 1);
    estimateStatT.Tr = zeros(nS, 1);
    
    
    %For caculating the Tr with the correct Sign and Value from EstimateR
    estimateValueR = abs(targetR).*sign(realTargetR);
    for i = 1:nS
        estimateValueStatS.Tr(i) = getTr(int2S(i,:), estimateValueR, refR);
        estimateValueStatR.Tr(i) = getTr(int2R(i,:), estimateValueR, refR);
        estimateValueStatT.Tr(i) = getTr(int2T(i,:), estimateValueR, refR);
    end
    estimateValueStatS.Tr = estimateValueStatS.Tr/sqrt(Len*(Len-1)/2);
    estimateValueStatR.Tr = estimateValueStatR.Tr/sqrt(Len*(Len-1)/2);
    estimateValueStatT.Tr = estimateValueStatT.Tr/sqrt(Len*(Len-1)/2);

    sortestimateValueStatT = sort(estimateValueStatT.Tr);
    estimateValueAbove95S = sum(estimateValueStatS.Tr>sortestimateValueStatT(int16(nS*0.95)));
    estimateValueAbove99S = sum(estimateValueStatS.Tr>sortestimateValueStatT(int16(nS*0.99)));
    
    %For caculating the Tr with EstimateR
    for i = 1:nS
        estimateStatS.Tr(i) = getTr(int2S(i,:), targetR, refR);
        estimateStatR.Tr(i) = getTr(int2R(i,:), targetR, refR);
        estimateStatT.Tr(i) = getTr(int2T(i,:), targetR, refR);
    end
    estimateStatS.Tr = estimateStatS.Tr/sqrt(Len*(Len-1)/2);
    estimateStatR.Tr = estimateStatR.Tr/sqrt(Len*(Len-1)/2);
    estimateStatT.Tr = estimateStatT.Tr/sqrt(Len*(Len-1)/2);

    sortestimateStatT = sort(estimateStatT.Tr);
    estimateAbove95S = sum(estimateStatS.Tr>sortestimateStatT(int16(nS*0.95)));
    estimateAbove99S = sum(estimateStatS.Tr>sortestimateStatT(int16(nS*0.99)));


    %For caculating the Tr with correct Sign and correct R value
    for i = 1:nS
        correctStatS.Tr(i) = getTr(int2S(i,:), realTargetR, refR);
        correctStatR.Tr(i) = getTr(int2R(i,:), realTargetR, refR);
        correctStatT.Tr(i) = getTr(int2T(i,:), realTargetR, refR);
    end
    correctStatS.Tr = correctStatS.Tr/sqrt(Len*(Len-1)/2);
    correctStatR.Tr = correctStatR.Tr/sqrt(Len*(Len-1)/2);
    correctStatT.Tr = correctStatT.Tr/sqrt(Len*(Len-1)/2);
    
    sortCorrectStatT = sort(correctStatT.Tr);
    correctAbove95S = sum(correctStatS.Tr>sortCorrectStatT(int16(nS*0.95)));
    correctAbove99S = sum(correctStatS.Tr>sortCorrectStatT(int16(nS*0.99)));
    
    %For caculating the Tr befor recombination (with Sign from Ref and
    %estimate value
    for i = 1:nS
        preStatS.Tr(i) = getTr(int2S(i,:), preTargetR, refR);
        preStatR.Tr(i) = getTr(int2R(i,:), preTargetR, refR);
        preStatT.Tr(i) = getTr(int2T(i,:), preTargetR, refR);
    end
    preStatS.Tr = preStatS.Tr/sqrt(Len*(Len-1)/2);
    preStatR.Tr = preStatR.Tr/sqrt(Len*(Len-1)/2);
    preStatT.Tr = preStatT.Tr/sqrt(Len*(Len-1)/2);
    
    sortPreStatT = sort(preStatT.Tr);
    preAbove95S = sum(preStatS.Tr>sortPreStatT(int16(nS*0.95)));
    preAbove99S = sum(preStatS.Tr>sortPreStatT(int16(nS*0.99)));
    
    %finalTargetR = preTargetR;
    
    %for blocks larger than 3
    blocks = [blocks blocks(:,2) - blocks(:,1) + 1];
    [m n] = size(blocks);
    %add learning marker
    blocks = [blocks, zeros(m, 1)];
    %blocks = blocks(blocks(:,3)>=3, :);
    
    [m n] = size(blocks);
    finalResult = zeros(m, m);
    
    f = fopen('result.txt', 'w');
    
    %Save all the sign matrix
    bufferMatrix = zeros(Len, Len, trials);
    
    %% inner block learning

    currentSeq = innerBlockDriver(caseSeq4, refSeq4, blocks);
    %currentSeq = refSeq4;
 
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
                block1(1,4) = 1;
                %currentSeq = shuffleNewBlock(currentSeq, block2);
                parfor t = 1:trials

				%profile on;
                    [finalSeq finalR finalSignRate finalQual blockMask] = newHBRecombo(targetR, caseSeq4, refSeq4, newCurrentSeq, block1, block2, alleleMapping, 0);

                    % p = profile('info');
                    % profsave(p,'profile_results');

                    bufferMatrix(:,:,t) = sign(finalR.*blockMask);
                    blockRate(t,:) = [finalSignRate finalQual];
                    if finalQual == 0 && finalSignRate ~= 1.0
                        pause
                    end
                end
            else
                block2(1,4) = 1;
                %currentSeq = shuffleNewBlock(currentSeq, block1);
                  parfor t = 1:trials

					%profile on;
                    [finalSeq finalR finalSignRate finalQual blockMask] = newHBRecombo(targetR, caseSeq4, refSeq4, newCurrentSeq, block2, block1, alleleMapping, 0);

                    % p = profile('info');
                    % profsave(p,'profile_results');

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
            
            fprintf(1, 'optimal r dist = %f, signRate = %f between block %d and block %d\n', blockRate(minIdx,2),blockRate(minIdx,1), i, j);
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
        postStatT.Tr(i) = getTr(int2T(i,:), finalTargetR, refR);
    end
    
    postStatS.Tr = postStatS.Tr/sqrt(Len*(Len-1)/2);
    postStatR.Tr = postStatR.Tr/sqrt(Len*(Len-1)/2);
    postStatT.Tr = postStatT.Tr/sqrt(Len*(Len-1)/2);
    
    sortPostStatT = sort(postStatT.Tr);
    postAbove95S = sum(postStatS.Tr>sortPostStatT(int16(nS*0.95)));
    postAbove99S = sum(postStatS.Tr>sortPostStatT(int16(nS*0.99)));
    
    %% Plot all the results
    
    plotScatterNew(caseSeq4, refSeq4, testSeq4, preTargetR, refR, 'WithSignsfromREF', iBigRepeat);
    plotScatterNew(caseSeq4, refSeq4, testSeq4, finalTargetR, refR, 'AfterSignRecovery', iBigRepeat);
    plotScatterNew(caseSeq4, refSeq4, testSeq4, targetR, refR, 'WithEstimateR', iBigRepeat);
	plotScatterNew(caseSeq4, refSeq4, testSeq4, realTargetR, refR, 'WithCalcR', iBigRepeat);
	
	[preSignRate test1 test2] = SignRate(realTargetR, refR);
	[postSignRate test1 test2] = SignRate(realTargetR, finalTargetR);
	[estiSignRate test1 test2] = SignRate(realTargetR, targetR);
	[postToEstiSignRate test1 test2] = SignRate(finalTargetR, targetR);
    [preToEstiSignRate test1 test2] = SignRate(refR, targetR);
    
    
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
    averagePreEstiSign = averagePreEstiSign + preToEstiSignRate;
	
	averagePre95 = averagePre95 + preAbove95S;
	averagePost95 = averagePost95 + postAbove95S;
	averageCorrect95 = averageCorrect95 + correctAbove95S;
	averageEsti95 = averageEsti95 + estimateValueAbove95S;
	averageHomer95 = averageHomer95 + homerAbove95S;
    averageRealEsti95 = averageRealEsti95 + estimateAbove95S;
	
	averagePre99 = averagePre99 + preAbove99S;
	averagePost99 = averagePost99 + postAbove99S;
	averageCorrect99 = averageCorrect99 + correctAbove99S;
	averageEsti99 = averageEsti99 + estimateValueAbove99S;
	averageHomer99 = averageHomer99 + homerAbove99S;
    averageRealEsti99 = averageRealEsti99 + estimateAbove99S;
    
    fprintf (1, 'PreSignRate = %f\n ', preSignRate);
    fprintf (1, 'PostSignRate = %f\n ', postSignRate);
	fprintf (1, 'EstiSignRate = %f\n ', estiSignRate);
	fprintf (1, 'PostToEstiSignRate = %f\n ', postToEstiSignRate);
	fprintf (1, 'PreToEstiSignRate = %f\n ', preToEstiSignRate);
	
    fprintf (1, '**************************************\n');
    fprintf (1, 'PreAbove95s             %d\n ', preAbove95S);
    fprintf (1, 'PostAbove95s            %d\n ', postAbove95S);
    fprintf (1, 'CorrectAbove95s         %d\n ', correctAbove95S);
    fprintf (1, 'CorrectSignEstimateRAbove95s            %d\n ', estimateValueAbove95S);
    fprintf (1, 'EstimateRAbove95S       %d\n ', estimateAbove95S);
    fprintf (1, 'HomerAbove95s           %d\n ', homerAbove95S);
	
	fprintf (1, '**************************************\n');
    fprintf (1, 'PreAbove99s             %d\n ', preAbove99S);
    fprintf (1, 'PostAbove99s            %d\n ', postAbove99S);
    fprintf (1, 'CorrectAbove99s         %d\n ', correctAbove99S);
    fprintf (1, 'CorrectSignEstimateRAbove99s            %d\n ', estimateValueAbove99S);
    fprintf (1, 'EstimateRAbove99S       %d\n ', estimateAbove99S);
    fprintf (1, 'HomerAbove99s           %d\n ', homerAbove99S);
%     fprintf (1, '**************************************\n');
%     fprintf (1, 'PreCaseTr>0         %d\n ', sum(preStatS.Tr>0));
%     fprintf (1, 'PostCaseTr>0        %d\n ', sum(postStatS.Tr>0));
%     fprintf (1, 'PreReferenceTr<0    %d\n ', sum(preStatR.Tr<0));
%     fprintf (1, 'PostReferenceTr<0   %d\n ', sum(postStatR.Tr<0));
%     fprintf (1, 'CorrectCaseTr>0     %d\n ', sum(correctStatS.Tr>0));
%     fprintf (1, 'CorrectReference<0  %d\n ', sum(correctStatR.Tr<0));

fprintf (1, '**************************************\n');
fprintf (1, 'Average Test results for %d bigRepeat and %d trials\n' , iBigRepeat, trials);
fprintf (1, '**************************************\n');
fprintf (1, 'AveragePreSignRate = %f\n ', averagePreSign/iBigRepeat);
fprintf (1, 'AveragePostSignRate = %f\n ', averagePostSign/iBigRepeat);
fprintf (1, 'AverageEstiSignRate = %f\n ', averageEstiSign/iBigRepeat);
fprintf (1, 'AveragePreToEstiSignRate = %f\n ', averagePreEstiSign/iBigRepeat);
fprintf (1, 'AveragePostToEstiSignRate = %f\n ', averagePostEstiSign/iBigRepeat);
fprintf (1, '**************************************\n');
fprintf (1, 'AveragePreAbove95s             %f\n ', averagePre95/iBigRepeat);
fprintf (1, 'AveragePostAbove95s            %f\n ', averagePost95/iBigRepeat);
fprintf (1, 'AverageCorrectAbove95s         %f\n ', averageCorrect95/iBigRepeat);
fprintf (1, 'AverageRealAbove95s         %f\n ', averageRealEsti95/iBigRepeat);
fprintf (1, 'AverageCorrectSignEstimateRAbove95s            %f\n ', averageEsti95/iBigRepeat);
fprintf (1, 'AverageHomerAbove95s           %f\n ', averageHomer95/iBigRepeat);
fprintf (1, '**************************************\n');
fprintf (1, 'AveragePreAbove99s             %f\n ', averagePre99/iBigRepeat);
fprintf (1, 'AveragePostAbove99s            %f\n ', averagePost99/iBigRepeat);
fprintf (1, 'AverageCorrectAbove99s         %f\n ', averageCorrect99/iBigRepeat);
fprintf (1, 'AverageRealAbove99s         %f\n ', averageRealEsti99/iBigRepeat);
fprintf (1, 'AverageCorrectSignEstimateRAbove99s            %f\n ', averageEsti99/iBigRepeat);
fprintf (1, 'AverageHomerAbove99s           %f\n ', averageHomer99/iBigRepeat);
end

fprintf (1, '**************************************\n');
fprintf (1, 'Average Test results for %d bigRepeat and %d trials\n' , bigRepeat, trials);
fprintf (1, '**************************************\n');
fprintf (1, 'AveragePreSignRate = %f\n ', averagePreSign/bigRepeat);
fprintf (1, 'AveragePostSignRate = %f\n ', averagePostSign/bigRepeat);
fprintf (1, 'AverageEstiSignRate = %f\n ', averageEstiSign/bigRepeat);
fprintf (1, 'AveragePreToEstiSignRate = %f\n ', averagePreEstiSign/bigRepeat);
fprintf (1, 'AveragePostToEstiSignRate = %f\n ', averagePostEstiSign/bigRepeat);
fprintf (1, '**************************************\n');
fprintf (1, 'AveragePreAbove95s             %f\n ', averagePre95/bigRepeat);
fprintf (1, 'AveragePostAbove95s            %f\n ', averagePost95/bigRepeat);
fprintf (1, 'AverageCorrectAbove95s         %f\n ', averageCorrect95/bigRepeat);
fprintf (1, 'AverageRealAbove95s         %f\n ', averageRealEsti95/bigRepeat);
fprintf (1, 'AverageCorrectSignEstimateRAbove95s            %f\n ', averageEsti95/bigRepeat);
fprintf (1, 'AverageHomerAbove95s           %f\n ', averageHomer95/bigRepeat);
fprintf (1, '**************************************\n');
fprintf (1, 'AveragePreAbove99s             %f\n ', averagePre99/bigRepeat);
fprintf (1, 'AveragePostAbove99s            %f\n ', averagePost99/bigRepeat);
fprintf (1, 'AverageCorrectAbove99s         %f\n ', averageCorrect99/bigRepeat);
fprintf (1, 'AverageRealAbove99s         %f\n ', averageRealEsti99/bigRepeat);
fprintf (1, 'AverageCorrectSignEstimateRAbove99s            %f\n ', averageEsti99/bigRepeat);
fprintf (1, 'AverageHomerAbove99s           %f\n ', averageHomer99/bigRepeat);

fclose(f);
diary off;

save;
save('finalResult.mat', 'finalResult');
save('caseSeq4.mat', 'caseSeq4');
save('refSeq4.mat', 'refSeq4');
save('finalTargetR.mat', 'finalTargetR');
