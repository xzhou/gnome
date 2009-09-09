function [finalCaseSeq, finalCaseR, signRate] = simHyplotypeRecombo(caseSeq4, refSeq4, blocks)
%this function use the hyplotype recombination technology to recover the
%sign of case R
    
    %get the major allele mapping skeme
    alleleMapping = getMajorAllele(refSeq4);
    
    %calculate R value
    caseR = calcR(caseSeq4, alleleMapping);
    refR = calcR(refSeq4, alleleMapping);
    
    %start matlabpool for parallel computing
    numOfCPU = 2;
    isOpen = matlabpool('size') > 0;
    if ~isOpen
        matlabpool(numOfCPU);
    end
    
    nBlocks = length(blocks(:1));
    
    %do recombination for each pair of blocks
    for i = 1:nBlocks-1
        for j = i:nBlocks
            aPairOfBlocks = blocks([i j],:);
            [targetR refSeq blockAlleleMapping] = getTarget(caseR, refSeq4, alleleMapping, aPairOfBlocks);
            
            result = cell(20, 3);
            parfor k = 1:20
                [finalSeq finalR finalSignRate] = hbOneBlockRecombination(targetR, refSeq, blockAlleleMapping, aPairOfBlocks);
                result(k,:) = {finalSeq; finalR; finalSignRate};
            end
            
            %select the best singRate
            signRates = [result{:,3}];
            [maxSignRate index] = max(signRates);
            bestSeq = result{index, 1};
            bestR = result{index,2};
            %TODO add this to the caseSequence
        end
    end
end