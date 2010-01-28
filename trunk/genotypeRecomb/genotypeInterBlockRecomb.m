function [result] = genotypeInterBlockRecomb(caseSeq, startingSeq, config)
    %recombination of genotype sequence
    blocks = config.blocks;
    [nBlocks tmp] = size(blocks);
    blocks = [blocks blocks(:,2) - blocks(:,1) + 1];
    blocks = [blocks, zeros(nBlocks, 1)];
    
    %==============config==================
    trials = config.trials;
    nInterBlockRecomb = config.nInterBlockRecomb;
    alpha = config.alpha;    %the weight of single allele frequency
    %==============end config==============
    
    [nIndividual nSnps] = size(startingSeq);

    bestR = estimateR(startingSeq);
    targetR = estimateR(caseSeq);
        
    learningResult = zeros(nSnps, nSnps, trials);
    for i = 1:(nBlocks-1)
        for j = (i+1):nBlocks    
            blockRate = zeros(nInterBlockRecomb, 2);
            block1 = blocks(i,:);
            block2 = blocks(j,:);
            %for each block, try many times
            for itr = 1:trials
                if block1(1,3) >= block2(1,3)
                    block1(1, 4) = 1;
                    for t = 1:nInterBlockRecomb
                        [finalSeq finalR finalSignRate finalQual blockMask] = pairGenotypeRecombo(targetR, caseSeq, startingSeq, block1, block2, config);
                        %merge the reuslt
                        learningResult(:, :, itr) = sign(finalR.*blockMask);
                        blockRate(itr, :) = [finalSignRate, finalQual];
                    end
                else
                    block2(1, 4) = 1;
                    for t = 1:nInterBlockRecomb
                        [finalSeq finalR finalSignRate finalQual blockMask] = pairGenotypeRecombo(targetR, caseSeq, startingSeq, block2, block1, config);
                        %merge the result
                        learningResult(:,:,itr) = sign(finalQual.*blockMask);
                        blockRate(itr, :) = [finalSignRate, finalQual];
                    end
                end
            end
            
            [minQ, minIdx] = min(blockRate(:,2));
            maxQ = min(blockRate(:,2));
            
            %assign the optimal sign to final TargetR
            bestR = bestR.*(learningResult(:,:,minIdx)==0) + abs(bestR).*(learningResult(:,:,minIdx)~=0).*learningResult(:,:,minIdx);
        end
    end
    
    result.finalTargetR = bestR;
end
