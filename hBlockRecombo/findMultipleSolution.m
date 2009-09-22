%% find multiple solutions for target R
function [finalR, finalRDiff, finalQuality] = findMultipleSolution(caseSeq4, refSeq4, blocks)
    debugMode = 1;
    %calculation
    alleleMapping = getMajorAllele(refSeq4);
    refR = calcR(refSeq4, alleleMapping);
    caseR = calcR(caseSeq4, alleleMapping);
    
    %add the size of the block to the blocks
    [nBlock nBlockCol] = size(blocks);
    blocks = [blocks blocks(:,2)-blocks(:,1) + 1 zeros(nBlock, 1)];
    
    firstBlock = getNextBlock(blocks);
    blocks(firstBlock, 4) = 1;
    
    %learn from current Sequence
    currentSeq = refSeq4;
    
    i = getNextBlock(blocks);

    newStart = blocks(i, 1);
    newEnd = blocks(i, 2);

    newBlock = blocks(i,:);
    
    solutions = zeros(10,1);
    
    startParallel(2);
    
    
    for i = 1:10
        [tempFinalSeq tempFinalR tempFinalSignRate finalQual] = newHBRecombo(caseR, caseSeq4, currentSeq, blocks, newBlock, alleleMapping);
        diff = caseSeq4 == tempFinalSeq;
        nDiff = sum(sum(diff));
       if finalQual == 0 && nDiff ~= 0
            printf(1, 'multiple solution found!!!\n');
            solutions(i) = 1;
        end
    end

end