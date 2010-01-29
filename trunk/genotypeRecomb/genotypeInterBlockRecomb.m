function [result] = genotypeInterBlockRecomb(caseSeq, startingSeq, config)
    %recombination of genotype sequence
    disp 'starting inter block recombination'
    blocks = config.blocks;
    [nBlocks tmp] = size(blocks);
    blocks = [blocks blocks(:,2) - blocks(:,1) + 1];
    blocks = [blocks, zeros(nBlocks, 1)];
    
    %==============config==================
    trials = config.trials;
    alpha = config.alpha;    %the weight of single allele frequency
    %==============end config==============
    
    [nIndividual nSnps] = size(startingSeq);

    bestR = estimateR(startingSeq);
    targetR = estimateR(caseSeq);
        
    learningResult = zeros(nSnps, nSnps, trials);
    for i = 1:(nBlocks-1)
        for j = (i+1):nBlocks    
            blockRate = zeros(trials, 2);
            block1 = blocks(i,:);
            block2 = blocks(j,:);
            %for each block, try many times
            if block1(1,3) >= block2(1,3)
                block1(1, 4) = 1;
                parfor t = 1:trials
                    [finalSeq finalR finalSignRate finalQual blockMask] = pairGenotypeRecombo(targetR, caseSeq, startingSeq, block1, block2, config);
                    %merge the reuslt
                    learningResult(:, :, t) = sign(finalR.*blockMask);
                    blockRate(t, :) = [finalSignRate, finalQual];
                end
            else
                block2(1, 4) = 1;
                parfor t = 1:trials
                    [finalSeq finalR finalSignRate finalQual blockMask] = pairGenotypeRecombo(targetR, caseSeq, startingSeq, block2, block1, config);
                    %merge the result
                    learningResult(:,:,t) = sign(finalQual.*blockMask);
                    blockRate(t, :) = [finalSignRate, finalQual];
                end
            end
            
            [minQ, minIdx] = min(blockRate(:,2));
            maxQ = min(blockRate(:,2));
            %assign the optimal sign to final TargetR
            bestR = bestR.*(learningResult(:,:,minIdx)==0) + abs(bestR).*(learningResult(:,:,minIdx)~=0).*learningResult(:,:,minIdx);
        end
    end
    
    result.finalTargetR = bestR;
    save('interblock');
end
