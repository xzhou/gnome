%this is the main function
function [] = WTCCC1_1500_Genotype_Recombo()
    addpath '~/research_linux/gnome/bioWorkspace/genomeprj/common';
    change_env()    %change the environment
    startParallel(); %start parallelel
    
    %>>>>>>>>>>>>>>> start configuration >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    wtccc1Conf.dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/TPED';
    wtccc1Conf.genotypeFile = 'Affx_gt_58C_Chiamo_07.tped.extract.inp.ped.fixed';
    wtccc1Conf.phaseFastaFile = 'Affx_gt_58C_Chiamo_07.tped.fasta';
    wtccc1Conf.blocks = [1, 24; 25, 65; 66, 81];       %manually define the strucuture
    wtccc1Conf.sampleSize = 250;    %the number of individuals for case or reference
    wtccc1Conf.logFileName = 'WTCCC1_1500.log';
    
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
    
    %% init r
    %[totalR pA counts] = estimateR(genotypeAll);
    
    %% for experiments
    [caseSeq, refSeq, caseID, refID] = randomSelectGenotype(genotypeAll, idInfo, wtccc1Conf.sampleSize);
    casePhaseSeq = getPhaseSeq(caseID, haplotypeSeq);
    refPhaseSeq = getPhaseSeq(refID, haplotypeSeq);
    
    
    
    %analysis seeverqever
    casePhaseIntSeq = seq2int(casePhaseSeq);
    refPahseIntSeq = seq2int(refPhaseSeq);
    [caseRef, summaryResult, caseInfo, refInfo] = analysisPhasedCaseRef(casePhaseIntSeq, refPahseIntSeq, wtccc1Conf.blocks);
    save('caseref.mat', 'caseSeq', 'refSeq', 'caseID', 'refID', 'casePhaseSeq', 'refPhaseSeq', 'caseRef', 'summaryResult');
        
    [m n] = size(caseSeq);
    fprintf(wtccc1Conf.logfid, 'sample size %d X %d\n', m, n);

    %% doing innerblock learning to approach the frequency
    fprintf(wtccc1Conf.logfid, 'start inner block learning');
    [randomCaseSeq] = gInnerSeqLearning(caseSeq, refSeq, wtccc1Conf);
    
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