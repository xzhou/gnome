%% do relaxed hyplotype block recombination for real data. 
%We will try to relax the consistency model for free recombination between
%any blocks to see if we can recover more signs;

%start computing slaves
startParallel(2);

blocks = [1 1; 2 2; 3 3;4 4; 5 5; 6 6; 7 7; 8 8; 9 9; 10 16; 
            18 21;22 22; 23 23; 24 24; 25 28; 29 29; 30 30; 31 31; 32 32; 33 33; 
            34 41; 41 41; 42 42; 43 43; 44 44; 45 45; 46 46; 47 47; 48 48; 49 49; 50 53; 
            54 56; 57 57; 58 58; 59 59; 60 60; 61 61; 62 62; 63 63; 64 64; 65 65; 66 66;
            67 67; 68 68; 69 69; 70 70;
            71 77];
       
%blocks = [1 15; 16 59; 60 77];            

blocks = [1 24; 25 45; 46 111; 112 174];      


%cd 'D:\IUBResearch\Projects\Bioinfor\data\88_77_CEU_YRI_DATA';
cd 'D:\IUBResearch\Projects\Bioinfor\data\HAPMAP';
%cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP';

[nBlock tmp] = size(blocks);

delete('hbrecombo.log');
diary hbrecombo.log;

%rawFastaData = fastaread('hapmap_chr7_80SNP_CEU_haplotype.fasta');
rawFastaData = fastaread('chr10_FGFR2_200kb_phased_CEU.fasta');

for iBigRepeat = 1:1
    fprintf(1, '\n*** trial %d ***\n', iBigRepeat);
    
    
%      for 77 SNPs randomly select case and reference
    [caseSeq4 refSeq4] = randomSelect(rawFastaData);
    
    
    
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
    homerAbove95S = sum(StatS.Tp>sortHomerStatR(int8(nS*0.95)));
    figure;
    hold on;
    plot(index1, StatS.Tp, '.r');
    plot(index2, StatR.Tp, '.g');
    legend({'case' 'ref'});
    plot(ones(2*nS).*sortHomerStatR(int8(nS*0.95)));
    xlabel('individual index');
    ylabel('T_r value');
    title('Homer Test');


    targetR = calcR(caseSeq4, alleleMapping);
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
    correctAbove95S = sum(correctStatS.Tr>sortCorrectStatR(int8(nS*0.95)));

    %For caculating the Tr befor recombination
    for i = 1:nS
        preStatS.Tr(i) = getTr(int2S(i,:), preTargetR, refR);
        preStatR.Tr(i) = getTr(int2R(i,:), preTargetR, refR);
    end
    preStatS.Tr = preStatS.Tr/sqrt(Len*(Len-1)/2); 
    preStatR.Tr = preStatR.Tr/sqrt(Len*(Len-1)/2);

    sortPreStatR = sort(preStatR.Tr);
    preAbove95S = sum(preStatS.Tr>sortPreStatR(int8(nS*0.95)));

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


    %initCaseSeq = currentSeq;
    %% Generate start point accoring to ref
    newCurrentSeq = zeros(nS, Len);
    currentBlockFreqInfo = cell(nBlock, 1);
    parfor i = 1:nBlock
        currentBlockFreqInfo{i,1} = getBlockFreq(currentSeq, blocks(i,:));
    end
    currentMatchedRef = blockCheck(currentBlockFreqInfo, refBlockFreqInfo, blocks);
    currentMatchedCase = blockCheck(currentBlockFreqInfo, caseBlockFreqInfo, blocks);

    for i = 1:nBlock      
        currentMatchedRef{i, 1} = currentMatchedRef{i,1}(1:end-2, :);
        for j = 1 : nS
            tempRefSeq4 = refSeq4(:, blocks(i,1):blocks(i,2));
            k = 1; 
            found = 0;
            while ((k<=length(currentMatchedRef{i,1}(:,1))&&(found == 0)))
                found = isequal(tempRefSeq4(j,:), currentMatchedRef{i,1}(k, 1:end-2));
                k = k+1;
            end
            if ((found == 1)&&(currentMatchedRef{i,1}(k-1, end-1)~=0))
                newCurrentSeq(j, blocks(i,1):blocks(i,2)) = currentMatchedRef{i,1}(k-1, 1:end-2);
                currentMatchedRef{i,1}(k-1, end-1) = currentMatchedRef{i,1}(k-1, end-1)-1;
                currentMatchedRef{i,1}(k-1, end) = currentMatchedRef{i,1}(k-1, end)-1;
            else
                disBetweenSeq = zeros(length(currentMatchedRef{i,1}), 1);
                for l = 1 : length(currentMatchedRef{i,1}(:,1))
                    if ((currentMatchedRef{i,1}(l,end-1))>(currentMatchedRef{i,1}(l, end)))
                        disBetweenSeq(l,1) = sum(tempRefSeq4(j,:) == currentMatchedRef{i,1}(l, 1:end-2));
                    end
                end
                [maxDisSum,maxDisIdx] = max(disBetweenSeq(:,1));
                newCurrentSeq(j, blocks(i,1):blocks(i,2)) = currentMatchedRef{i,1}(maxDisIdx, 1:end-2);
                currentMatchedRef{i,1}(maxDisIdx, end-1) = currentMatchedRef{i,1}(maxDisIdx, end-1)-1;
            end
        end
    end


    finalTargetR = calcR(currentSeq, alleleMapping);
    

    for i = 1:(m-1)
        for j = i+1:m
            blockRate = zeros(trials, 2);      
            block1 = blocks(i,:);
            block2 = blocks(j,:);
            if(block1(1,3) >= block2(1,3))
                block1(1,4) = 1;
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
                block2(1,4) = 1;
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
    postAbove95S = sum(postStatS.Tr>sortPostStatR(int8(nS*0.95)));

    plotScatter(caseSeq4, refSeq4, preTargetR, refR, 'With Signs from REF');
    plotScatter(caseSeq4, refSeq4, finalTargetR, refR, 'After Sign Recovery');
    plotScatter(caseSeq4, refSeq4, targetR, refR, 'With Correct Signs');
    preSignRate = sum(sum(sign(targetR)==sign(refR)))/Len/(Len-1);
    postSignRate = sum(sum(sign(targetR)==sign(finalTargetR)))/Len/(Len-1);
    
    
    fprintf (1, ' PreSignRate = %f\n ', preSignRate);
    fprintf (1, 'PostSignRate = %f\n ', postSignRate);

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
