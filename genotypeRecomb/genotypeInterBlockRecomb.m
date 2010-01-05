function [result] = genotypeInterBlockRecomb(caseSeq, startingSeq, blocks, alleleMapping, false)
    %recombination of genotype sequence
    [nBlocks tmp] = size(blocks);
    if tmp <= 2
      %extend the block to pre calculate the size
      blocks = [blocks blocks(:,2) - blocks(:,1)]
    end

    %==============config==================
    trials = 1000;
    nInterBlockRecomb = 1000;
    alpha = 0.01    %the weight of 
    %==============end config==============
    
    [nIndividual nSnps] = size(startingSeq);

    initR = estimateR(startingSeq, alleleMapping);

    result = zeros(nSnps, nSnps, trials);
    
    for itr = 1:trials
        for i = 1:(nBlocks-1)
            for j = (i+1):nBlocks
                block1 = blocks(i,:);
                block2 = blocks(j,:);
                if block1(1,3) >= block2(1,3)
                    for t = 1:nInterBlockRecomb
                        [finalSeq finalR finalSignRate finalQual blockMask] = pairGenotypeRecombo(targetR, caseSeq, startingSeq, block1, block2, alleleMapping, 0.01);
                        %merge the reuslt
                    end
                else
                    for t = 1:nInterBlockRecomb
                        [finalSeq finalR finalSignRate finalQual blockMask] = pairGenotypeRecombo(targetR, caseSeq, startingSeq, block2, block1, alleleMapping, 0.01);
                        %merge the result
                    end
                end
            end
        end
    end
    
    
  
end
