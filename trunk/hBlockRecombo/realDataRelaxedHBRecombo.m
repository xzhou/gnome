%% do relaxed hyplotype block recombination for real data. 
%We will try to relax the consistency model for free recombination between
%any blocks to see if we can recover more sings;

blocks = [1 1; 2 2; 3 3;4 4; 5 5; 6 6; 7 7; 8 8; 9 9; 10 16; 
            18 21;22 22; 23 23; 24 24; 25 28; 29 29; 30 30; 31 31; 32 32; 33 33; 
            34 41; 41 41; 42 42; 43 43; 44 44; 45 45; 46 46; 47 47; 48 48; 49 49; 50 53; 
            54 56; 57 57; 58 58; 59 59; 60 60; 61 61; 62 62; 63 63; 64 64; 65 65; 66 66;
            67 67; 68 68; 69 69; 70 70;
            71 77];

cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP';

delete('hbrecombo.log');
diary hbrecombo.log;

real80SNPSeq = readSeq4('hapmap_chr7_80SNP_CEU_haplotype.fasta');

[m n] = size(real80SNPSeq);

caseSeq4 = real80SNPSeq(1:m/2, :);
refSeq4 = real80SNPSeq((m/2+1):end, :);

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

f = fopen('result.txt', 'w');

for i = 1:(m-1)
    for j = i+1:m
        blockRate = zeros(trials, 2);            
        block1 = blocks(i,:);
        block2 = blocks(j,:);            
        if(block1(1,3) >= block2(1,3))
            block1(1,4) = 1;
            parfor t = 1:trials
                [finalSeq finalR finalSignRate finalQual] = newHBRecombo(targetR, caseSeq4, refSeq4, block1, block2, alleleMapping);
                blockRate(t,:) = [finalSignRate finalQual];
            end
        else
            block2(1,4) = 1;
            parfor t = 1:trials
                [finalSeq finalR finalSignRate finalQual] = newHBRecombo(targetR, caseSeq4, refSeq4, block2, block1, alleleMapping);
                blockRate(t,:) = [finalSignRate finalQual];
            end
        end

        [maxVal, maxIdx] = max(blockRate(:,1));
        maxQ = max(blockRate(:,2));
        
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

fclose(f);
diary off;

save('finalResult.mat', 'finalResult');
