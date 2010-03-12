%this is the main function
function [] = ALL_GENOTYPE_WTCCC1_1500_Genotype_Recombo()
%this function will start from a block pool and first sample each block and
%then connect them
    addpath '~/research_linux/gnome/bioWorkspace/genomeprj/common';
    change_env()    %change the environment
    startParallel(); %start parallelel
    
    %% ===== start configuration ===========
    wtccc1Conf.dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/TPED';
    wtccc1Conf.genotypeFile = 'Affx_gt_58C_Chiamo_07.tped.extract.inp.ped.fixed';
    wtccc1Conf.phaseFastaFile = 'Affx_gt_58C_Chiamo_07.tped.fasta';
    wtccc1Conf.blocks = [1, 24; 25, 65; 66, 81];       %manually define the strucuture
    wtccc1Conf.invariantSnps = [];
    wtccc1Conf.caseSize= 250;    %the number of individuals for case or reference
    wtccc1Conf.refSize = 1250;   %we use more sample in reference group
    wtccc1Conf.logFileName = 'RELAX_HFREQ_HBSTART_WTCCC1_1500_Genotype_Recombo.log';

    %output options
    wtccc1Conf.verbose = true;
    
    wtccc1Conf.highFreqT = 5;
    
    %inner block learning configuration
    wtccc1Conf.innerBlockExpT = 1.0e-6;
    wtccc1Conf.maxItr = 10000;
    wtccc1Conf.nRepeat = 10;
    
    %inter block learning configuration
    wtccc1Conf.trials = 10;
    wtccc1Conf.nInterBlockRecomb = 1;
    wtccc1Conf.alpha = 0.01;    %combination of r^2 diff and cx0 idff
    wtccc1Conf.smallFilter = 0;
    
    wtccc1Conf.maxIT = 10000;
    %======= end configuration ==================
    
    %===== initialization =============
    if wtccc1Conf.verbose
        wtccc1Conf.logfid = 1;
    else
        wtccc1Conf.logfid = fopen(wtccc1Conf.logFileName, 'w');
    end
    %=======< end initialization ===========>>
    
    %% goto data directory
    cd(wtccc1Conf.dataPath);
    [genotypeAll majorAllele, idInfo] = readPedFile(wtccc1Conf.genotypeFile);
    disp(['reading ', wtccc1Conf.genotypeFile, ' complete']);
    
    [haplotypeSeq] = fastaread(wtccc1Conf.phaseFastaFile);
    %allHapSeqNoID = removeHapID(haplotypeSeq);
    %allHapIntSeqNoID = seq2int(allHapSeqNoID);
    
    
    %% for experiments
    [caseSeq, refSeq, caseID, refID] = randomSelectGenotype(genotypeAll, idInfo, wtccc1Conf);
    
    %get phased sequence
    casePhaseSeq = getPhaseSeq(caseID, haplotypeSeq);
    refPhaseSeq = getPhaseSeq(refID, haplotypeSeq);

    %analysis sequences
    casePhaseIntSeqWithID = seq2int(casePhaseSeq);
    refPahseIntSeqWithID = seq2int(refPhaseSeq);
    casePhaseIntSeqNoID = getSeqMatrix(casePhaseIntSeqWithID);
    refPhaseIntSeqNoID = getSeqMatrix(refPahseIntSeqWithID);
    save('caseRefInit.mat');
    
    %analysis target estimiate R and R
    casePhaseR = calcR(casePhaseIntSeqNoID, majorAllele);
    caseEstR = estimateR(caseSeq);
    signRateMax = SignRate(casePhaseR, caseEstR);
    fprintf(wtccc1Conf.logfid, 'case est r sign rate (MAX sign rate) = %f', signRateMax);
    save('calcRComplete.mat');
    
    [~, blockCoverRate, caseFreqInfo, refFreqInfo] = analysisPhasedCaseRef(casePhaseIntSeqNoID, refPhaseIntSeqNoID, wtccc1Conf.blocks);
    printCoverRate(wtccc1Conf.logfid, blockCoverRate, 'random hyplotype');
    
    [m n] = size(caseSeq);
    fprintf(wtccc1Conf.logfid, 'sample size %d X %d\n', m, n);
    
    
    %% do data analysis
%     [sampledHapSeq] = sampleHapSeq(refPhaseIntSeqNoID, wtccc1Conf);
%     [sampledHapSeq] = blockSampleHapSeq(refPhaseIntSeqNoID, wtccc1Conf);
%     [~, sampleHapCoverRate, ~, ~] = analysisPhasedCaseRef(casePhaseIntSeqNoID, sampledHapSeq, wtccc1Conf.blocks);
%     printCoverRate(wtccc1Conf.logfid, sampleHapCoverRate, 'sample hap cover rate');
%     
%     [sampledGenoSeq] = hapSeq2GenoSeq(sampledHapSeq, majorAllele);
%     [~, sampleCoverRate, caseFreqInfo, sampleFreqInfo] = analysisPhasedCaseRef(caseSeq, sampledGenoSeq, wtccc1Conf.blocks);
%     printCoverRate(wtccc1Conf.logfid, sampleCoverRate, 'sample genotype cover rate');
%     
    randomRef = randomKRow(refSeq, wtccc1Conf.caseSize);
    [randomSummary, randomCoverRate, caseFreqInfo, randomFreqInfo] = analysisPhasedCaseRef(caseSeq, randomRef, wtccc1Conf.blocks);
    printCoverRate(wtccc1Conf.logfid, randomCoverRate, 'random seq cover rate');
    %plotCoverRate(caseFreqInfo, 'case');
%     
%     %total cover rate analysis
%     [totalSummary, totalCoverRate, ~, totalRefFreqInfo] = analysisPhasedCaseRef(caseSeq, refSeq, wtccc1Conf.blocks);
%     printCoverRate(wtccc1Conf.logfid, totalCoverRate, 'total seq cover rate');
%     plotCoverRate(totalRefFreqInfo, 'total freq info');

    %%shuffle case frequencies, small distubance
    save;
    [shuffleRef] = smallDisturbance(caseFreqInfo);
    [shuffleSummary, shuffleCoverRate, ~, shuffleFreqInfo] = analysisPhasedCaseRef(caseSeq, shuffleRef, wtccc1Conf.blocks);
    printCoverRate(wtccc1Conf.logfid, shuffleCoverRate, 'shuffle');
    plotCoverRate(shuffleFreqInfo, 'shuffle');
    plotCoverRate(caseFreqInfo, 'case');
    
    %% doing innerblock learning to approach the frequency
    fprintf(wtccc1Conf.logfid, 'start inner block learning');
    [randomCaseSeq] = gInnerSeqLearning(caseSeq, shuffleRef, wtccc1Conf);
    
    caseSeqAfterInnerBlockLearning = adjustStartPoint(randomCaseSeq, refSeq, wtccc1Conf.blocks);
    
    %% doing interblock learning starting from new sequence
    [result] = genotypeInterBlockRecomb(caseSeq, caseSeqAfterInnerBlockLearning, wtccc1Conf);

    %check the sign recover rate against phasing since we don't have the
    %real data
    %calculate target estimated R
    targetR = estimateR(caseSeq);
    initRefR = estimateR(refSeq);
     
    initSignRate = SignRate(targetR, initRefR);
    finalSignRate = SignRate(targetR, result.finalTargetR);

    fprintf(1, 'initSignRate %f, \tfinaSignRate, %f', initSignRate, finalSignRate);
    
end