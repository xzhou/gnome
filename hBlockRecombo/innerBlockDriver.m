function [seq] = innerBlockDriver(caseSeq4, refSeq4, blocks, verbose)    
%cd 'D:\IUBResearch\Projects\Bioinfor\data\88_77_CEU_YRI_DATA';
%     clear;
   % cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP';
    
    if nargin == 0
        rawFastaData = fastaread('hapmap_chr7_80SNP_CEU_haplotype.fasta');
        [caseSeq4 refSeq4] = randomSelect(rawFastaData);    
        blocks = [1 15; 16 55; 60 77];
    end
    
    if nargin == 3
        verbose = true;
    end
 
    [nBlock tmp] = size(blocks);

    nS = length(caseSeq4);
    Len = length(caseSeq4(1,:));

    caseBlockFreqInfo = cell(nBlock, 1);

    parfor i = 1:nBlock
        caseBlockFreqInfo{i,1} = getBlockFreq(caseSeq4, blocks(i,:));
    end

    refBlockFreqInfo = cell(nBlock, 1);
    parfor i = 1:nBlock
        refBlockFreqInfo{i,1} = getBlockFreq(refSeq4, blocks(i,:));
    end

    matchedCase = blockCheck(caseBlockFreqInfo, refBlockFreqInfo, blocks);
    refMatchedCase = blockCheck(refBlockFreqInfo, caseBlockFreqInfo, blocks);

    startParallel();

    nRepeat = 10;
    result = cell(nBlock,nRepeat);

    for i=1:nBlock
        block = blocks(i,:);
        a = block(1, 1);
        b = block(1, 2);
        refBlock = refMatchedCase{i,1}(1:end-2,1:end-2);
        refFreq = refMatchedCase{i,1}(1:end-2, end-1);
        refCaseFreq = refMatchedCase{i,1}(1:end-2, end);
        caseBlock = matchedCase{i,1}(1:end-2, 1:end-2);
        caseFreq = matchedCase{i,1}(1:end-2, end-1);
        caseRefFreq = matchedCase{i,1}(1:end-1, end);
        refSeq = refSeq4(:,a:b);
        caseSeq = caseSeq4(:,a:b);
        alleleMapping = getMajorAllele(refSeq);
        if verbose
            fprintf(1, '\n************ learning block %d ***********\n', i);
        end
        for k = 1:nRepeat
            profile on;
            [aResult] = innerBlockLearning(caseBlock, caseFreq, refBlock, refFreq, refCaseFreq, alleleMapping);
                                    profile viewer;
                    p = profile('info');
                    profsave(p,'profile_results')
            %[bResult] = innerBlockLearning(caseBlock, caseFreq, refBlock, refFreq, refCaseFreq, alleleMapping,1);
            try
                result{i, k} = aResult;
            catch exception
                fprintf(1, 'error\n');
                rethrow(exception);
            end
            if verbose
                fprintf(1, 'block = %d repeat = %d a = %f initSR = %f, finalSR = %f\n',i, k, aResult.fDistance, aResult.initSignRate, aResult.finalSignRate);
            end
        end
    end
    %save('result.mat', 'result');
    
    %% reconstruct the learned block we use the
    seq = InnerBlockHelp.recoverCaseSeq(result);
    save('innterseq.mat', 'seq');
    save;
end