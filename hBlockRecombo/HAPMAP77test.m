%in this test, we are going to use CEU hapmap data and the generate a
%markov chain model and generate multiple haplotype sequence for
%experiments

cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP'

%for each pair of blocks, we repeat the experiments for repeatTimes to get
%a better result.
repeatTimes = 100;

%ceuSequence = readSeq4('chr10_FGFR2_200kb_phased_CEU.fasta');
ceuSeqStruct = fastaread('chr10_FGFR2_200kb_phased_CEU.fasta');
ceuSeq4 = readSeq4('chr10_FGFR2_200kb_phased_CEU.fasta');
yriSeqStruct = fastaread('chr10_FGFR2_200kb_phased_yri.fasta');
yriSeq4 = readSeq4('chr10_FGFR2_200kb_phased_yri.fasta');

%get the major allele mapping
alleleMapping = getMajorAllele(ceuSeq4);
ceuSeq01 = (ceuSeq4 == repmat(alleleMapping, length(ceuSeq4), 1)) + 0;
yriSeq01 = (yriSeq4 == repmat(alleleMapping, length(yriSeq4), 1)) + 0;

[mCEU nCEU] = size(ceuSeq4);
[mYRI nYRI] = size(ceuSeq4);

if 0
    if loadModel == 0
        %generate the markov chain model
        iMCmodelCEU = iMCmodelBuild(ceuSeqStruct, 0);
        iMCmodelYRI = iMCmodelBuild(yriSeqStruct, 0);
        save('iMCmodelCEU.mat', 'iMCmodelCEU');
        save('iMCmodelYRI.mat', 'iMCmodelYRI');
    else
        load('iMCmodelCEU.mat');
        load('iMCmodelYRI.mat');
    end
end

%directly use the read data for the learning
%%show LD similarity
ceuR = calcR(ceuSeq4, alleleMapping);
yriR = calcR(yriSeq4, alleleMapping);

% plotWithSignSimple(abs(ceuR), abs(yriR), 0, 0);
% xlabel('yri');
% ylabel('ceu');
% saveas(gcf, 'ceuR_yriR.pdf');

%%cut the graph into 3 parts for easy process
absCeuR = ceuR;
absYriR = yriR;
plotWithSignSimple(absCeuR(1:80, 1:80), absYriR(1:80,1:80), 0, 0);
xlabel('yri 1:80');
ylabel('ceu 1:80');
saveas(gcf, 'ceuR_yriR 1_80.pdf');
% 
% plotWithSignSimple(absCeuR(50:120, 50:120), absYriR(50:120, 50:120), 0, 0);
% xlabel('yri 50:120');
% ylabel('ceu 50:120');
% saveas(gcf, 'ceuR_yriR 50_120.pdf');
% 
% plotWithSignSimple(absCeuR(90:end, 90:end), absYriR(90:end, 90:end), 0, 0);
% xlabel('yri 90:end');
% ylabel('ceu 90:end');
% saveas(gcf, 'ceuR_yriR 90_end.pdf');



%start matlab pool for parallel compuation
isOpen = matlabpool('size') > 0;
if ~isOpen
    matlabpool 2;
end

%identifiy blocks
blocks = [1 10; 12 23; 26 44];

%%learn first block
% [targetR refSeq blockAlleleMapping] = getTarget(ceuR, yriSeq4, alleleMapping, blocks(1:2,:));
% parfor i = 1:20
%     [finalSeq finalR finalSignRate] = hbOneBlockRecombination(targetR, refSeq, blockAlleleMapping, blocks(1:2,:));
% end

%do the experiments 20 times
blocks = [1 23; 25 46];
[targetR refSeq blockAlleleMapping] = getTarget(ceuR, yriSeq4, alleleMapping, blocks(1:2,:));
parfor i = 1:20
    [finalSeq finalR finalSignRate] = hbOneBlockRecombination(targetR, refSeq, blockAlleleMapping, blocks(1:2,:));
end

