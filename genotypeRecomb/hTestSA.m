%% do relaxed hyplotype block recombination for real data. 
%We will try to relax the consistency model for free recombination between
%any blocks to see if we can recover more signs;

try
    addpath '/home/xzhou/research_linux/gnome/workspace/hBlockRecombo'
    addpath '/home/xzhou/research_linux/gnome/workspace/MyDrTest'
catch exception
    %do nothing
end

%start computing slaves

startParallel();

blocks = [1, 24; 25, 65; 66, 81]; 

% blocks = [1 24; 25 45; 46 111; 112 174];      


%cd 'D:\IUBResearch\Projects\Bioinfor\data\88_77_CEU_YRI_DATA';
%for different platform

try
    cd 'D:\IUBResearch\Projects\Bioinfor\data\1500';
    disp 'WINDOWS'
catch e
    disp 'LINUX'
    cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP';
end


[nBlock tmp] = size(blocks);

delete('hbrecombo.log');
diary hbrecombo.log;

%rawFastaData = fastaread('hapmap_chr7_80SNP_CEU_haplotype.fasta');
rawFastaData = fastaread('newAffx.fasta');

for iBigRepeat = 1:1
    fprintf(1, '\n*** trial %d ***\n', iBigRepeat);
    
%      for 77 SNPs randomly select case and reference
    [caseSeq4 refSeq4] = randomSelect(rawFastaData, 250, 0);
    
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


    %targetRealR = calcR(caseSeq4, alleleMapping);
    targetR = calcR(caseSeq4, alleleMapping);
    realSingleAlleleFreq = GnomeCalculator.getSingleAlleleFreq(caseSeq4, alleleMapping);
    %round to precision 1e-4
    
    %targetGenotypeSeq = haplotype2genotype(caseSeq4, alleleMapping);
    %[targetR pA counts] = estimateR(targetGenotypeSeq);
    targetR = fix(targetR.*10000)./10000;
    
    realTargetR = targetR;
    
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

    trials = 16;

    f = fopen('result.txt', 'w');

    %Save all the sign matrix
    bufferMatrix = zeros(Len, Len, trials);

    %% inner block learning
    %currentSeq1 = innerBlockDriver(caseSeq4, refSeq4, blocks);
    
    type innerBlockMutate.m
    type innerBlockFitness.m
    currentSeq = zeros(size(refSeq4));
    m=nBlock;
    for i = 1:m
        %fitnessfcn = @(x) innerBlockFitness(x,caseSeq(:,1:15));
        tempCaseR = calcR(caseSeq4(:,blocks(i,1):blocks(i,2)), alleleMapping(1, blocks(i,1):blocks(i,2)));
        tempCaseRs = tempCaseR.*tempCaseR;
        tempCaseAlleleFreq = GnomeCalculator.getSingleAlleleFreq(caseSeq4(:,blocks(i,1):blocks(i,2)), alleleMapping(1, blocks(i,1):blocks(i,2)));
        fitnessfcn = @(x) innerBlockFitness(x,tempCaseRs, tempCaseAlleleFreq, alleleMapping(1, blocks(i,1):blocks(i,2)));
        options = saoptimset('DataType', 'custom', 'AnnealingFcn', @innerBlockMutate, ...
            'StallIterLimit',10000, 'ReannealInterval', 10000, 'TolFun', 1e-70);
        % Finally, we call simulated annealing with our problem information.
        [genotypeBlock,funval, exitflag, output] = simulannealbnd(fitnessfcn, refSeq4(:,blocks(i,1):blocks(i,2)), [], [], options);
        %genotypeBlock = simulannealbnd(fitnessfcn, refSeq, [], [], options);
        
        currentSeq(:,blocks(i,1):blocks(i,2)) = genotypeBlock;
    end   
    
    r1 = calcR(currentSeq);
    %r11 = calcR(currentSeq1);
    r2 = calcR(caseSeq4);
    r3 = calcR(refSeq4);
    diff1 = sum(sum(abs(r1.*r1-r2.*r2)))/2;
    diff2 = sum(sum(abs(r2.*r2-r3.*r3)))/2;
%     diff3 = sum(sum(abs(r11.*r11-r2.*r2)))/2;
      
    currentBlockFreqInfo = cell(nBlock, 1);
    parfor i = 1:nBlock
        currentBlockFreqInfo{i,1} = getBlockFreq(currentSeq, blocks(i,:));
    end
    currentMatchedRef = blockCheck(currentBlockFreqInfo, refBlockFreqInfo, blocks);
    currentMatchedCase = blockCheck(currentBlockFreqInfo, caseBlockFreqInfo, blocks);
    
%     current1BlockFreqInfo = cell(nBlock, 1);
%     parfor i = 1:nBlock
%         current1BlockFreqInfo{i,1} = getBlockFreq(currentSeq1, blocks(i,:));
%     end
%     current1MatchedRef = blockCheck(current1BlockFreqInfo, refBlockFreqInfo, blocks);
%     current1MatchedCase = blockCheck(current1BlockFreqInfo, caseBlockFreqInfo, blocks);
      
    newCurrentSeq = zeros(nS, Len);
    newCurrentSeq = adjustStartPoint(currentSeq, refSeq4, blocks);
    
    finalTargetR = calcR(newCurrentSeq, alleleMapping);
    finalTargetR = abs(targetR).*sign(finalTargetR);
    
    %% Start the inter block recombination
    
%    newCurrentSeq = refSeq4;
    
    type interBlockMutate.m
    type interBlockFitness.m
    
    signMatrix = zeros(size(caseSeq4));
%   seqAfterInnerLearning = haplotype2genotype(currentSeq);
    signMatrix = sign(finalTargetR);
    for i= 1:(m-1)
        for j = i+1:m
            block1 = blocks(i,:);
            block2 = blocks(j,:);
            if(block1(1,3)>= block2(1,3))
                block1(1,4) = 1;
                fitnessfcn = @(x) interBlockFitness(x,caseSeq4,block1,block2);
                interBlockMutateNew = @(optimValues, problem) interBlockMutate(optimValues, problem, block2);
                options = saoptimset('DataType', 'custom', 'AnnealingFcn', interBlockMutateNew, ...
                    'StallIterLimit',10000, 'ReannealInterval', 10000, 'TolFun', 1e-70);
                genotypeBlock = simulannealbnd(fitnessfcn, newCurrentSeq, [], [], options);
                
            else
                block2(1,4) = 1;
                fitnessfcn = @(x) interBlockFitness(x,caseSeq4,block2,block1);
                interBlockMutateNew = @(optimValues, problem) interBlockMutate(optimValues, problem, block1);
                options = saoptimset('DataType', 'custom', 'AnnealingFcn', interBlockMutateNew, ...
                    'StallIterLimit',10000, 'ReannealInterval', 10000, 'TolFun', 1e-70);
                genotypeBlock = simulannealbnd(fitnessfcn, newCurrentSeq, [], [], options);
            end
            blockMask = getBlockMaskForEval(finalTargetR, block1, block2);
            signMatrix = signMatrix.*(blockMask==0) + sign(calcR(genotypeBlock)).*blockMask;  
        end
        %Get the sign from the seq after learning at the end
 
    end
    
    finalTargetR = signMatrix.*finalTargetR;
    
    

    %% For caculating the Tr after recombination
    for i = 1:nS
        postStatS.Tr(i) = getTr(int2S(i,:), finalTargetR, refR);
        postStatR.Tr(i) = getTr(int2R(i,:), finalTargetR, refR);
    end

    postStatS.Tr = postStatS.Tr/sqrt(Len*(Len-1)/2);
    postStatR.Tr = postStatR.Tr/sqrt(Len*(Len-1)/2);

    sortPostStatR = sort(postStatR.Tr);
    postAbove95S = sum(postStatS.Tr>sortPostStatR(int16(nS*0.95)));

%     plotScatter(caseSeq4, refSeq4, preTargetR, refR, 'With Signs from REF');
%     plotScatter(caseSeq4, refSeq4, finalTargetR, refR, 'After Sign Recovery');
%     plotScatter(caseSeq4, refSeq4, targetR, refR, 'With Correct Signs');
%     preSignRate = sum(sum(sign(targetR)==sign(refR)))/Len/(Len-1);
%     postSignRate = sum(sum(sign(targetR)==sign(finalTargetR)))/Len/(Len-1);

	[preSignRate test1 test2] = SignRate(realTargetR, refR);
	[postSignRate test1 test2] = SignRate(realTargetR, finalTargetR);
	[estiSignRate test1 test2] = SignRate(realTargetR, targetR);
	[postToEstiSignRate test1 test2] = SignRate(finalTargetR, targetR);
    
    
    fprintf (1, ' PreSignRate = %f\n ', preSignRate);
    fprintf (1, 'PostSignRate = %f\n ', postSignRate);

%     fprintf (1, 'PreCaseTr>0         %d\n ', sum(preStatS.Tr>0));
%     fprintf (1, 'PostCaseTr>0        %d\n ', sum(postStatS.Tr>0));
%     fprintf (1, 'PreReferenceTr<0    %d\n ', sum(preStatR.Tr<0));
%     fprintf (1, 'PostReferenceTr<0   %d\n ', sum(postStatR.Tr<0));
%     fprintf (1, 'CorrectCaseTr>0     %d\n ', sum(correctStatS.Tr>0));
%     fprintf (1, 'CorrectReference<0  %d\n ', sum(correctStatR.Tr<0));
end


fclose(f);
diary off;

save;
save('finalResult.mat', 'finalResult');
save('caseSeq4.mat', 'caseSeq4');
save('refSeq4.mat', 'refSeq4');
save('finalTargetR.mat', 'finalTargetR');
