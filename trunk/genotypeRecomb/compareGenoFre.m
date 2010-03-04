function [] = compareGenoFre(caseSeq, sampledGenoSeq, wtccc1Conf)

    blocks = wtccc1Conf.blocks;
    [nBlock tmp] = size(blocks);
    
    caseBlockFreqInfo = cell(nBlock, 1);

    parfor i = 1:nBlock
    caseBlockFreqInfo{i,1} = getBlockFreq(caseSeq, blocks(i,:));
    end

    caseSampleBlockFreqInfo = cell(nBlock, 1);
    
    parfor i = 1:nBlock
    caseSampleBlockFreqInfo{i,1} = getBlockFreq(sampledGenoSeq, blocks(i,:));
    end

    caseMatchedCaseSample = blockCheck(caseBlockFreqInfo, caseSampleBlockFreqInfo, blocks);
    caseSampleMatchedCase = blockCheck(caseSampleBlockFreqInfo, caseBlockFreqInfo, blocks);
    
end