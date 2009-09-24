%We will try to relax the consistency model for free recombination between
%any blocks to see if we can recover more sings;

blocks = [1 2 3 4 5 6 7 9  12 18 22 23 24 30 31 32 34 42 46 50 51 53 54 58 61 65 68 69 71 72 73 74; 
          1 2 3 4 5 6 8 11 17 21 22 23 29 30 31 33 41 45 49 50 52 53 57 60 64 67 68 70 71 72 73 77]';

%blocks = [34 41; 42 45];
      
cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77'

refFileName = 'SIM_100x77_ctl.fasta';
sampleFileName = 'SIM_100x77_smp.fasta';


caseSeq4 = readSeq4(sampleFileName);
refSeq4 = readSeq4(refFileName);

alleleMapping = getMajorAllele(refSeq4);

targetR = calcR(caseSeq4, alleleMapping);
refR = calcR(refSeq4, alleleMapping);

%for blocks larger than 3
blocks = [blocks blocks(:,2) - blocks(:,1) + 1];
[m n] = size(blocks);
%add learning marker
blocks = [blocks, zeros(m, 1)];
blocks = blocks(blocks(:,3)>=3, :);

startParallel(2);
[m n] = size(blocks);
finalResult = zeros(m, m);

trials = 10;

f = open('result.txt');

for i = 1:m-1
    for j = i+1:m
        blockRate = zeros(trials, 1);            
        block1 = blocks(i,:);
        block2 = blocks(j,:);            
        if(block1(1,3) >= block2(1,3))
            block1(1,4) = 1;
            parfor t = 1:trials
                [finalSeq finalR finalSignRate finalQual] = newHBRecombo(targetR, caseSeq4, refSeq4, block1, block2, alleleMapping);
                blockRate(t) = finalSignRate;
            end
        else
            block2(1,4) = 1;
            parfor t = 1:trials
                [finalSeq finalR finalSignRate finalQual] = newHBRecombo(targetR, caseSeq4, refSeq4, block2, block1, alleleMapping);
                blockRate(t) = finalSignRate;
            end
        end
        finalResult(i,j) = max(blockRate);
        fprintf(f, '[%d-%d]x[%d-%d]\t = %f\n', block1(1,1), block1(1,2), block2(1,1), block2(1,2), finalResult(i,j));
    end
end

save('finalResult.mat', 'finalResult');
