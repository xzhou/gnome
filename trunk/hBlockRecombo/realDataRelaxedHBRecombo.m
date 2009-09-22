%% do relaxed hyplotype block recombination for real data. 
%We will try to relax the consistency model for free recombination between
%any blocks to see if we can recover more sings;

blocks = [10 17; 17 21; 25 28; 34 41; 50 53; 54 56; 71 77];

      
cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP';

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

for i = 1:m-1
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
        else
            fprintf(f, '[%d-%d]x[%d-%d]\t = %f\t NO\n', block1(1,1), block1(1,2), block2(1,1), block2(1,2), finalResult(i,j));
        end
        
      
    end
end

fclose(f);

save('finalResult.mat', 'finalResult');
