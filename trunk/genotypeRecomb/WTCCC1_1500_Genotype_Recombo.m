%this is the main function
function [] = WTCCC1_1500_Genotype_Recombo()
    change_env()    %change the environment
    startParallel(2); %start parallelel
    
    %======== start configuration =========================================
    dataPath = '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/TPED';
    genotypeFile = 'Affx_gt_58C_Chiamo_07.tped.extract.inp.ped';
    blocks = [1, 24; 25, 65; 66, 81];       %manually define the strucuture
    sampleSize = 250;
    %======== end configuration ===========================================
    
    %goto data directory
    cd(dataPath);
    [genotypeAll majorAllele] = readPedFile(genotypeFile);
    disp(['reading ', genotypeFile, ' complete']);
    
    %% init r
    [totalR pA counts] = estimateR(genotypeAll);
    
    %% for experiments
    [caseSeq, refSeq] = randomSelectGenotype(genotypeAll);
    
%     %calculate target estimated R
%     parfor pi = 0:1
%       if pi == 0
%         targetR = estimateR(caseSeq);
%       elseif pi == 1
%         initRefR = estimateR(refSeq);
%       end
%     end
    
    %% doing innerblock learning to approach the frequency
    disp(['start inner block learning']);
    [randomCaseSeq] = gInnerSeqLearning(caseSeq, refSeq, blocks, false);
    
    caseSeqAfterInnerBlockLearning = adjustStartPoint(randomCaseSeq, refSeq, blocks);
    
    %% doing interblock learning starting from new sequence
    [result] = genotypeInterBlockRecomb(caseSeq, caseSeqAfterInnerBlockLearning, blocks, false);
    
    %check the sign recover rate
end