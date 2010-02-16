%MAIN
function [] = SEQ_SEL_HBSTART_WTCCC1_1500_Genotype_Recombo()
% this function first select a reference genotype sequence that with
% closest rs difference and then adjust the genotype frequency and then
% reconstruct the genoteyp sequence
    addpath '~/research_linux/gnome/bioWorkspace/genomeprj/common';
    change_env()    %change the environment
    startParallel(); %start parallelel
    
    %>>>>>>>>>>>>>>> start configuration >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    wtccc1Conf.dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/TPED';
    wtccc1Conf.genotypeFile = 'Affx_gt_58C_Chiamo_07.tped.extract.inp.ped.fixed';
    wtccc1Conf.phaseFastaFile = 'Affx_gt_58C_Chiamo_07.tped.fasta';
    wtccc1Conf.blocks = [1, 24; 25, 65; 66, 81];       %manually define the strucuture
    wtccc1Conf.caseSize= 250;    %the number of individuals for case or reference
    wtccc1Conf.refSize = 1250;   %we use more sample in reference group
    wtccc1Conf.logFileName = 'WTCCC1_1500.log';

    %hyplotype pre selection configuration
    wtccc1Conf.alg = 1;
    wtccc1Conf.maxItr = 0;
    
    %output options
    wtccc1Conf.verbose = true;
    
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
    %<<<<<<<<<<<<<<< end configuration <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    %>>>>>>>>>>>>>>> initialization >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if wtccc1Conf.verbose
        wtccc1Conf.logfid = 1;
    else
        wtccc1Conf.logfid = fopen(wtccc1Conf.logFileName, 'w');
    end
    %<<<<<<<<<<<<<<< end initialization >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    
    %goto data directory
    cd(wtccc1Conf.dataPath);
    [genotypeAll majorAllele, idInfo] = readPedFile(wtccc1Conf.genotypeFile);
    disp(['reading ', wtccc1Conf.genotypeFile, ' complete']);
    
    [haplotypeSeq] = fastaread(wtccc1Conf.phaseFastaFile);

    %% for experiments
    [caseSeq, refSeq, caseID, refID] = randomSelectGenotype(genotypeAll, idInfo, wtccc1Conf);
    
    %get phased sequence
    casePhaseSeq = getPhaseSeq(caseID, haplotypeSeq);
    refPhaseSeq = getPhaseSeq(refID, haplotypeSeq);
    
    save('caseref.mat', 'caseSeq', 'refSeq', 'caseID', 'refID', 'casePhaseSeq', 'refPhaseSeq');
    
    %analysis sequences
    casePhaseIntSeqWithID = seq2int(casePhaseSeq);
    refPahseIntSeqWithID = seq2int(refPhaseSeq);
    casePhaseIntSeqNoID = getSeqMatrix(casePhaseIntSeqWithID);
    refPhaseIntSeqNoID = getSeqMatrix(refPahseIntSeqWithID);
    
    [hBlockSummary, blockCoverRate, caseFreqInfo, refFreqInfo] = analysisPhasedCaseRef(casePhaseIntSeqNoID, refPhaseIntSeqNoID, wtccc1Conf.blocks);
    
    [m n] = size(caseSeq);
    fprintf(wtccc1Conf.logfid, 'sample size %d X %d\n', m, n);
    
    [sampGenotype] = selectGenotype(refPhaseIntSeqNoID, wtccc1Conf);
    
    estTargetR = estimateR(caseSeq);
    estTargetRs = estTargetRs.*estTargetRs
    
    [refHapSeq] = sampleHapSeq(refPhaseIntSeqNoID, wtccc1Conf);
    sampledGenoSeq = getSmallDistanceSeqs(refHapSeq, m, estTargetRs, majorAllele);
    
    %[sampledHapSeq] = blockSampleHapSeq(refPhaseIntSeqNoID, wtccc1Conf);
    %[sampledGenoSeq] = hapSeq2GenoSeq(sampledHapSeq, majorAllele);
    
    %% doing innerblock learning to approach the frequency
    fprintf(wtccc1Conf.logfid, 'start inner block learning');
    [randomCaseSeq] = gInnerSeqLearning(caseSeq, sampledGenoSeq, wtccc1Conf);
    
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

    fprintf(1, 'initSignRate %f, \tfinalSignRate, %f', initSignRate, finalSignRate);
    
end