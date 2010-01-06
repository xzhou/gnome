function [result] = genotypeInterBlockRecomb(caseSeq, startingSeq, blocks, false)
    %recombination of genotype sequence
    [nBlocks tmp] = size(blocks);
    if tmp <= 2
      %extend the block to pre calculate the size
      blocks = [blocks blocks(:,2) - blocks(:,1)]
    end

    %==============config==================
    trials = 1000;
    nInterBlockRecomb = 1000;
    alpha = 0.01    %the weight of single allele frequency
    %==============end config==============
    
    [nIndividual nSnps] = size(startingSeq);

    initR = estimateR(startingSeq, alleleMapping);

    learningResult = zeros(nSnps, nSnps, trials);
    for i = 1:(nBlocks-1)
        for j = (i+1):nBlocks    
            blockRate = zeros(nInterBlockRecomb, 2);
            block1 = blocks(i,:);
            block2 = blocks(j,:);
            for itr = 1:trials
                if block1(1,3) >= block2(1,3)
                    for t = 1:nInterBlockRecomb
                        [finalSeq finalR finalSignRate finalQual blockMask] = pairGenotypeRecombo(targetR, caseSeq, startingSeq, block1, block2, 0.01);
                        %merge the reuslt
                        learningResult(:, :, itr) = sign(finalR.*blockMask);
                        blockRate(t, :) = [finalSignRate, finalQual];
                    end
                else
                    for t = 1:nInterBlockRecomb
                        [finalSeq finalR finalSignRate finalQual blockMask] = pairGenotypeRecombo(targetR, caseSeq, startingSeq, block2, block1, 0.01);
                        %merge the result
                        learningResult(:,:,itr) = sign(finalQual.*blockMask);
                        blockRate(t, :) = [finalSignRate, finalQual];
                    end
                end
            end
            
            [minQ, minIdx] = min(blockRate(:,2));
            maxQ = min(blockRate(:,2);
            
            finalTargetR = finalTargetR.*(bufferMatrix(:,:,minIdx)==0) + abs(finalTargetR).*(bufferMatrix(:,:,minIdx)~=0).*bufferMatrix(:,:,minIdx);
            finalResult(i,j) = maxVal;
        end
    end

end
