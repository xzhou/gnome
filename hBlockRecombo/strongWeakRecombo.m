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
        
blocks = [1 15; 16 55; 60 77];

[nBlock tmp] = size(blocks);

%cd 'D:\IUBResearch\Projects\Bioinfor\data\88_77_CEU_YRI_DATA';
cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP';

delete('hbrecombo.log');
diary hbrecombo.log;

rawFastaData = fastaread('hapmap_chr7_80SNP_CEU_haplotype.fasta');

[caseSeq4 refSeq4] = randomSelect(rawFastaData);
nS = length(caseSeq4);
Len = length(caseSeq4(1,:));

caseBlockFreqInfo = cell(nBlock, 1);

parfor i = 1:nBlock
	caseBlockFreqInfo{i,1} = getBlockFreq(caseSeq4, blocks(i,:));
end


%index for plotting
index1 = [1: nS];
index2 = [nS+1: nS*2];

alleleMapping = getMajorAllele(refSeq4);

int2S = (caseSeq4 == repmat(alleleMapping,nS,1)) + 0;
int2R = (refSeq4 == repmat(alleleMapping,nS,1)) + 0;

targetR = calcR(caseSeq4, alleleMapping);
refR = calcR(refSeq4, alleleMapping);

preTargetR = abs(targetR).*sign(refR);

preStatS.Tr = zeros(nS, 1);
preStatR.Tr = zeros(nS, 1);

postStatS.Tr = zeros(nS, 1);
postStatR.Tr = zeros(nS, 1);

correctStatS.Tr = zeros(nS, 1);
correctStatR.Tr = zeros(nS, 1);

h = plotScatter(caseSeq4, refSeq4, targetR, refR);


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

finalTargetR = preTargetR;



%for blocks larger than 3
blocks = [blocks blocks(:,2) - blocks(:,1) + 1];
[m n] = size(blocks);
%add learning marker
blocks = [blocks, zeros(m, 1)];
%blocks = blocks(blocks(:,3)>=3, :);


[m n] = size(blocks);
finalResult = zeros(m, m);

trials = 10;

f = fopen('result.txt', 'w');

%Save all the sign matrix
bufferMatrix = zeros(Len, Len, trials);

currentSeq = caseSeq4;


for i = 1:(m-1)
    for j = i+1:m
        blockRate = zeros(trials, 2);            
        block1 = blocks(i,:);
        block2 = blocks(j,:);
        if(block1(1,3) >= block2(1,3))
            block1(1,4) = 1;
            currentSeq = shuffleNewBlock(currentSeq, block2);
            parfor t = 1:trials
                [finalSeq finalR finalSignRate finalQual blockMask] = newHBRecombo(targetR, caseSeq4, currentSeq, block1, block2, alleleMapping, 0.01);
                bufferMatrix(:,:,t) = sign(finalR.*blockMask);
                blockRate(t,:) = [finalSignRate finalQual];
                if finalQual == 0 && finalSignRate ~= 1.0
                    pause
                end
            end
        else
            block2(1,4) = 1;
            currentSeq = caseSeq4;
            currentSeq = shuffleNewBlock(currentSeq, block1);
            parfor t = 1:trials
                [finalSeq finalR finalSignRate finalQual blockMask] = newHBRecombo(targetR, caseSeq4, currentSeq, block2, block1, alleleMapping, 0.01);
                bufferMatrix(:,:,t) = sign(finalR.*blockMask);
                blockRate(t,:) = [finalSignRate finalQual];
                
                if finalQual == 0 && finalSignRate ~= 1.0
                    pause
                end
            end
        end

        [maxVal, maxIdx] = max(blockRate(:,1));
        maxQ = min(blockRate(:,2));
        
        %Apply signs to these cross blocks
        finalTargetR = finalTargetR.*(bufferMatrix(:,:,maxIdx)==0) + abs(finalTargetR).*(bufferMatrix(:,:,maxIdx)~=0).*bufferMatrix(:,:,maxIdx);
        
        finalResult(i,j) = maxVal;
        
        if maxQ == blockRate(maxIdx, 2)
            fprintf(f, '[%d-%d]x[%d-%d]\t = %f\t YES\n', block1(1,1), block1(1,2), block2(1,1), block2(1,2), finalResult(i,j));
            fprintf(1, '[%d-%d]x[%d-%d]\t = %f\t YES\n', block1(1,1), block1(1,2), block2(1,1), block2(1,2), finalResult(i,j));
        else
            fprintf(f, '[%d-%d]x[%d-%d]\t = %f\t NO\n', block1(1,1), block1(1,2), block2(1,1), block2(1,2), finalResult(i,j));
            fprintf(1, '[%d-%d]x[%d-%d]\t = %f\t NO\n', block1(1,1), block1(1,2), block2(1,1), block2(1,2), finalResult(i,j));
        end
    end
end

%For caculating the Tr after recombination
for i = 1:nS
    postStatS.Tr(i) = getTr(int2S(i,:), finalTargetR, refR);
    postStatR.Tr(i) = getTr(int2R(i,:), finalTargetR, refR);
end

postStatS.Tr = postStatS.Tr/sqrt(Len*(Len-1)/2);

postStatR.Tr = postStatR.Tr/sqrt(Len*(Len-1)/2);

sortPostStatR = sort(postStatR.Tr);
postAbove95S = sum(postStatS.Tr>sortPostStatR(int8(nS*0.95)));




plotResult(index1, index2, preStatS.Tr, preStatR.Tr, postStatS.Tr, postStatR.Tr);
figure;
hold on;
plot(index1, correctStatS.Tr, '.r');
plot(index2, correctStatR.Tr, '.g');
preSignRate = sum(sum(sign(targetR)==sign(refR)))/77/77;
postSignRate = sum(sum(sign(targetR)==sign(finalTargetR)))/77/77;

fprintf (1, ' PreSignRate = %f\n ', preSignRate);
fprintf (1, 'PostSignRate = %f\n ', postSignRate);

fprintf (1, 'PreCaseTr>0         %d\n ', sum(preStatS.Tr>0));
fprintf (1, 'PostCaseTr>0        %d\n ', sum(postStatS.Tr>0));
fprintf (1, 'PreReferenceTr<0    %d\n ', sum(preStatR.Tr<0));
fprintf (1, 'PostReferenceTr<0   %d\n ', sum(postStatR.Tr<0));
fprintf (1, 'CorrectCaseTr>0     %d\n ', sum(correctStatS.Tr>0));
fprintf (1, 'CorrectReference<0  %d\n ', sum(correctStatR.Tr<0));


fclose(f);
diary off;

save('finalResult.mat', 'finalResult');
