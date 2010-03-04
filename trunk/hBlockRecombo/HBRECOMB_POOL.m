%% do relaxed hyplotype block recombination for real data. 
%We will try to relax the consistency model for free recombination between
%any blocks to see if we can recover more signs;

%change_env();
%start computing slaves

startParallel(2);

blocks = [1, 24; 25, 65; 66, 81]; 

try
    cd 'D:\IUBResearch\Projects\Bioinfor\data\1500';
    disp 'WINDOWS'
catch e
    try 
        cd '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/TPED';
        disp 'LINUX'
    catch e
    end
end

[nBlock tmp] = size(blocks);

delete('hbrecombo.log');
diary hbrecombo.log;

rawFastaData = fastaread('newAffx.fasta');
%rawFastaData = fastaread('Affx_gt_58C_Chiamo_07.tped.fasta');

for iBigRepeat = 1:1
    fprintf(1, '\n*** trial %d ***\n', iBigRepeat);
    
    
% Randomly select case and reference
    [caseSeq4 refSeq4pool] = randomSelect(rawFastaData, 250, 1);

    nS = length(caseSeq4(:,1));
    Len = length(caseSeq4(1,:));
    nR = length(refSeq4pool(:,1));



    %% test the block structure
    caseBlockFreqInfo = cell(nBlock, 1);

    parfor i = 1:nBlock
        caseBlockFreqInfo{i,1} = getBlockFreq(caseSeq4, blocks(i,:));
    end

    refPoolBlockFreqInfo = cell(nBlock, 1);
    parfor i = 1:nBlock
        refPoolBlockFreqInfo{i,1} = getBlockFreq(refSeq4pool, blocks(i,:));
    end
    
    caseMatchedRefPool = blockCheck(caseBlockFreqInfo, refPoolBlockFreqInfo, blocks);
    refPoolMatchedCase = blockCheck(refPoolBlockFreqInfo, caseBlockFreqInfo, blocks);
    

    for i = 1:nBlock
        [tempRows tempCols] = size(caseBlockFreqInfo{i,1});
        fprintf (1, '*********Block%d********* \n ', i);
        fprintf (1, 'Case    Hyplotype Blocks:%d \n ', tempRows);
        [tempRows tempCols] = size(refPoolBlockFreqInfo{i,1});
        fprintf (1, 'Ref     Hyplotype Blocks:%d \n ', tempRows);
        fprintf (1, 'Shared  Hyplotype Blocks:%d \n ', sum(caseMatchedRefPool{i,1}(:,end)~=0));
    end
    
    
    save('refBlockFreq.mat', 'refPoolBlockFreqInfo');
    save('caseBlockFreq.mat', 'caseBlockFreqInfo');
    save('caseMatchedRefPool.mat', 'caseMatchedRefPool');
    save('refPoolMatchedCase.mat', 'refPoolMatchedCase');
    
    %% Construct the reference using the pool
    refSeq4 = zeros(nS, Len);
    for i = 1:nBlock
        [tempRows tempCols] = size(refPoolBlockFreqInfo{i,1});
        refSeq4(1:tempRows,blocks(i,1):blocks(i,2)) = refPoolBlockFreqInfo{i,1}(:,1:end-1);
        refPoolBlockFreqInfo{i,1} = sortrows(refPoolBlockFreqInfo{i, 1}, -(blocks(i,2)-blocks(i,1)+2));
        j = 1;
        while (nS-tempRows)>round(refPoolBlockFreqInfo{i,1}(j,(blocks(i,2)-blocks(i,1)+2))/(nR/nS))
            for k = 1: round(refPoolBlockFreqInfo{i,1}(j,(blocks(i,2)-blocks(i,1)+2))/(nR/nS))
                tempRows = tempRows +1;
                refSeq4(tempRows,blocks(i,1):blocks(i,2)) = refPoolBlockFreqInfo{i,1}(j,1:(blocks(i,2)-blocks(i,1)+1));
            end
            j = j+1;
        end
        while (nS-tempRows) > 0
            tempRows = tempRows+1;
            refSeq4(tempRows,blocks(i,1):blocks(i,2)) = refPoolBlockFreqInfo{i,1}(j,1:(blocks(i,2)-blocks(i,1)+1));
        end
    end
    

    refBlockFreqInfo = cell(nBlock, 1);
    parfor i = 1:nBlock
        refBlockFreqInfo{i,1} = getBlockFreq(refSeq4, blocks(i,:));
    end
    
    matchedCase = blockCheck(caseBlockFreqInfo, refBlockFreqInfo, blocks);
    refMatchedCase = blockCheck(refBlockFreqInfo, caseBlockFreqInfo, blocks);


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

    targetR = calcR(caseSeq4, alleleMapping);
    realSingleAlleleFreq = GnomeCalculator.getSingleAlleleFreq(caseSeq4, alleleMapping);
    %round to precision 1e-4
    targetR = fix(targetR.*10000)./10000;
    
    refR = calcR(refSeq4, alleleMapping);

    preTargetR = abs(targetR).*sign(refR);

    preStatS.Tr = zeros(nS, 1);
    preStatR.Tr = zeros(nS, 1);

    postStatS.Tr = zeros(nS, 1);
    postStatR.Tr = zeros(nS, 1);

    correctStatS.Tr = zeros(nS, 1);
    correctStatR.Tr = zeros(nS, 1);


    %For caculating the Tr with correct Sign
    for i = 1:nS
        correctStatS.Tr(i) = getTr(int2S(i,:), targetR, refR);
        correctStatR.Tr(i) = getTr(int2R(i,:), targetR, refR);
    end
    correctStatS.Tr = correctStatS.Tr/sqrt(Len*(Len-1)/2);
    correctStatR.Tr = correctStatR.Tr/sqrt(Len*(Len-1)/2);


    sortCorrectStatR = sort(correctStatR.Tr);
    correctAbove95S = sum(correctStatS.Tr>sortCorrectStatR(int16(nS*0.95)));

    %For caculating the Tr befor recombination
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

    trials = 15;

    f = fopen('result.txt', 'w');

    %Save all the sign matrix
    bufferMatrix = zeros(Len, Len, trials);

    %% inner block learning
    currentSeq = innerBlockDriver(caseSeq4, refSeq4, blocks);
    
    currentBlockFreqInfo = cell(nBlock, 1);

    parfor i = 1:nBlock
        currentBlockFreqInfo{i,1} = getBlockFreq(currentSeq, blocks(i,:));
    end
    
    caseMatchedCurrent = blockCheck(caseBlockFreqInfo, currentBlockFreqInfo, blocks);
    currentMatchedCase = blockCheck(currentBlockFreqInfo, caseBlockFreqInfo, blocks);
    
    
    for i = 1:nBlock
        [tempRows tempCols] = size(caseBlockFreqInfo{i,1});
        fprintf (1, '*********Block%d********* \n ', i);
        fprintf (1, 'Case    Hyplotype Blocks:%d \n ', tempRows);
        [tempRows tempCols] = size(currentBlockFreqInfo{i,1});
        fprintf (1, 'Current Hyplotype Blocks:%d \n ', tempRows);
        fprintf (1, 'Shared  Hyplotype Blocks:%d \n ', sum(currentMatchedCase{i,1}(:,end)~=0));
    end

    newCurrentSeq = adjustStartPoint(currentSeq, refSeq4(1:nS,:), blocks);

    finalTargetR = calcR(newCurrentSeq, alleleMapping);
    
	finalTargetR = abs(targetR).*sign(finalTargetR);
    
    signMatrix = zeros(size(caseSeq4));
    signMatrix = sign(calcR(currentSeq));
    for i = 1:(m-1)
        for j = i+1:m
            blockRate = zeros(trials, 2);      
            block1 = blocks(i,:);
            block2 = blocks(j,:);
            if(block1(1,3) >= block2(1,3))
                block1(1,4) = 1;%??
                %currentSeq = shuffleNewBlock(currentSeq, block2);
                parfor t = 1:trials
                    [finalSeq finalR finalSignRate finalQual blockMask] = newHBRecombo(targetR, caseSeq4, newCurrentSeq, block1, block2, alleleMapping, 0.01);
                    bufferMatrix(:,:,t) = sign(finalR.*blockMask);
                    blockRate(t,:) = [finalSignRate finalQual];
                    if finalQual == 0 && finalSignRate ~= 1.0
                        pause
                    end
                end
            else
                block2(1,4) = 1;%??
                %currentSeq = shuffleNewBlock(currentSeq, block1);
                parfor t = 1:trials
                    [finalSeq finalR finalSignRate finalQual blockMask] = newHBRecombo(targetR, caseSeq4, newCurrentSeq, block2, block1, alleleMapping, 0.01);
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
%     plotScatter(caseSeq4, refSeq4, finalTargetR, refR,
%     'AfterSignRecovery', iBigRepeat);
    plotScatter(caseSeq4, refSeq4, targetR, refR, 'WithCorrectSigns', iBigRepeat);
    preSignRate = sum(sum(sign(targetR)==sign(refR)))/Len/(Len-1);
    postSignRate = sum(sum(sign(targetR)==sign(finalTargetR)))/Len/(Len-1);
    
    %preSignRate = sum(sum(sign(targetRealR)==sign(refR)))/Len/(Len-1);
    %postSignRate = sum(sum(sign(targetRealR)==sign(finalTargetR)))/Len/(Len-1);
    
    fprintf (1, ' PreSignRate = %f\n ', preSignRate);
    fprintf (1, 'PostSignRate = %f\n ', postSignRate);
    fprintf (1, '**************************************\n');
    fprintf (1, 'PreAbove95s             %d\n ', preAbove95S);
    fprintf (1, 'PostAbove95s            %d\n ', postAbove95S);
    fprintf (1, 'CorrectAbove95s         %d\n ', correctAbove95S);
    fprintf (1, 'HomerAbove95s           %d\n ', homerAbove95S);
    fprintf (1, '**************************************\n');
    fprintf (1, 'PreCaseTr>0         %d\n ', sum(preStatS.Tr>0));
    fprintf (1, 'PostCaseTr>0        %d\n ', sum(postStatS.Tr>0));
    fprintf (1, 'PreReferenceTr<0    %d\n ', sum(preStatR.Tr<0));
    fprintf (1, 'PostReferenceTr<0   %d\n ', sum(postStatR.Tr<0));
    fprintf (1, 'CorrectCaseTr>0     %d\n ', sum(correctStatS.Tr>0));
    fprintf (1, 'CorrectReference<0  %d\n ', sum(correctStatR.Tr<0));
end


fclose(f);
diary off;

save;
save('finalResult.mat', 'finalResult');
save('caseSeq4.mat', 'caseSeq4');
save('refSeq4.mat', 'refSeq4');
save('finalTargetR.mat', 'finalTargetR');
