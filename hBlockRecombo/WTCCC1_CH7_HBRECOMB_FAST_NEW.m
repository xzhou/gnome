function [] = WTCCC1_CH7_HBRECOMB_FAST_NEW()
%% do relaxed hyplotype block recombination for real data.
%We will try to relax the consistency model for free recombination between
%any blocks to see if we can recover more signs;

%OS = change_env();
%start computing slaves
startParallel();

addpath '~/GProject/code/hBlockRecombo/';
addpath '~/GProject/code/MyDrTest/';
addpath '~/GProject/code/genotypeRecomb/';
addpath '~/GProject/code/common/';

%*******************config********************%
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

wtccc1Conf.phaseFastaFile = 'newAffx.fasta';
wtccc1Conf.blocks = [1, 24; 25, 65; 66, 81];       %manually define the strucuture
wtccc1Conf.sampleSize = 250;    %the number of individuals for case or reference
wtccc1Conf.logFileName = 'WTCCC1_1500_HB_C00CONSTRAINT.log';

wtccc1Conf.nRep = 20;
wtccc1Conf.verbose = 1;

%*******************end config****************%

%*******************init**********************%
cd(wtccc1Conf.dataPath);
if wtccc1Conf.verbose
    wtccc1Conf.logfid = 1;
else
    wtccc1Conf.logfid = fopen(wtccc1Conf.logFileName, 'w');
end
%*******************end init******************%

[nBlock ~] = size(wtccc1Conf.blocks);
rawFastaData = fastaread(wtccc1Conf.phaseFastaFile);

[hap01Seq, alleleMapping] = encodeRawFastaSeq(rawFastaData);


for iBigRepeat = 1:wtccc1Conf.nRep
    fprintf(wtccc1Conf.logfid, '\n*** trial %d ***\n', iBigRepeat);
    
    [caseSeq4 refSeq4] = randomSelect(rawFastaData, 250);
    
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
        caseBlockFreqInfo{i,1} = getBlockFreq(caseSeq4, wtccc1Conf.blocks(i,:));
    end
    
    refBlockFreqInfo = cell(nBlock, 1);
    parfor i = 1:nBlock
        refBlockFreqInfo{i,1} = getBlockFreq(refSeq4, wtccc1Conf.blocks(i,:));
    end
    
    matchedCase = blockCheck(caseBlockFreqInfo, refBlockFreqInfo, wtccc1Conf.blocks);
    refMatchedCase = blockCheck(refBlockFreqInfo, caseBlockFreqInfo, wtccc1Conf.blocks);
    
    save('refBlockFreq.mat', 'refBlockFreqInfo');
    save('caseBlockFreq.mat', 'caseBlockFreqInfo');
    save('caseMatch.mat', 'matchedCase');
    save('refMatchedCase.mat', 'refMatchedCase');
    
    
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
    figure;
    hold on;
    plot(index1, StatS.Tp, '.r');
    plot(index2, StatR.Tp, '.g');
    legend({'case' 'ref'});
    plot(ones(2*nS).*sortHomerStatR(int16(nS*0.95)));
    xlabel('individual index');
    ylabel('T_r value');
    title('Homer Test');
    
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
    wtccc1Conf.blocks = [wtccc1Conf.blocks wtccc1Conf.blocks(:,2) - wtccc1Conf.blocks(:,1) + 1];
    [m n] = size(wtccc1Conf.blocks);
    %add learning marker
    wtccc1Conf.blocks = [wtccc1Conf.blocks, zeros(m, 1)];
    %blocks = blocks(blocks(:,3)>=3, :);
    
    
    [m n] = size(blocks);
    finalResult = zeros(m, m);
    
    trials = 15;
    
    f = fopen('result.txt', 'w');
    
    %Save all the sign matrix
    bufferMatrix = zeros(Len, Len, trials);
    
    %% inner block learning
    currentSeq = innerBlockDriver(caseSeq4, refSeq4, wtccc1Conf.blocks);
    
    newCurrentSeq = adjustStartPoint(currentSeq, refSeq4, wtccc1Conf.blocks);
    
    finalTargetR = calcR(currentSeq, alleleMapping);
    
    signMatrix = zeros(size(caseSeq4));
    signMatrix = sign(calcR(currentSeq));
    for i = 1:(m-1)
        for j = i+1:m
            blockRate = zeros(trials, 2);
            block1 = wtccc1Conf.blocks(i,:);
            block2 = wtccc1Conf.blocks(j,:);
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
    
    plotScatter(caseSeq4, refSeq4, preTargetR, refR, 'With Signs from REF');
    plotScatter(caseSeq4, refSeq4, finalTargetR, refR, 'After Sign Recovery');
    plotScatter(caseSeq4, refSeq4, targetR, refR, 'With Correct Signs');
    preSignRate = sum(sum(sign(targetR)==sign(refR)))/Len/(Len-1);
    postSignRate = sum(sum(sign(targetR)==sign(finalTargetR)))/Len/(Len-1);
    
    %preSignRate = sum(sum(sign(targetRealR)==sign(refR)))/Len/(Len-1);
    %postSignRate = sum(sum(sign(targetRealR)==sign(finalTargetR)))/Len/(Len-1);
    
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
end