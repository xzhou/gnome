%cd 'D:\IUBResearch\Projects\Bioinfor\data\88_77_CEU_YRI_DATA';
cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP';
blocks = [1 15; 16 55; 60 77];
[nBlock tmp] = size(blocks);

rawFastaData = fastaread('hapmap_chr7_80SNP_CEU_haplotype.fasta');
[caseSeq4 refSeq4] = randomSelect(rawFastaData);
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

startParallel(2);

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
    fprintf(1, 'block %d \n***********\n', i);
    for k = 1:nRepeat
        [aResult] = innerBlockLearning(caseBlock, caseFreq, refBlock, refFreq, refCaseFreq, alleleMapping);
        [bResult] = innerBlockLearning(caseBlock, caseFreq, refBlock, refFreq, refCaseFreq, alleleMapping,1);
        result{i, k} = aResult;
        fprintf(1, 'block = %d repeat = %d \t a = %f,\t b = %f\n',i, k, aResult.fDistance, bResult.fDistance);
    end
end

save('result.mat', 'result');

