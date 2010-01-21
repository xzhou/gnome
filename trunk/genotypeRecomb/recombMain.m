%this is the main function
function [] = recombMain()
    
    change_env()

    startParallel(2);
    
    
    %======== start configuration ==============
    
    
    dataPath = '/home/xzhou/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/TPED';
    genotypeFile = 'Affx_gt_58C_Chiamo_07.tped.extract.inp.ped';
    %manually defien the strucuture
    blocks = [1, 24; 25, 65; 66; 81];
    
    %======== end configuration   ==============
    
    all = readSeq4(fastaFile);
    realR = calcR(all);
    
    genotypeAll = genotypeHelpFuncs.readGenotypeFromFasta(fastaFile);
    
    %% init r
    [totalR pA counts] = estimateR(genotypeAll);
    
    %% for experiments
    [caseSeq, refSeq] = randomSelectGenotype(genotypeAll);
    
    %calculate target estimated R
%     parfor pi = 0:1
%       if pi == 0
%         targetR = estimateR(caseSeq);
%       elseif pi == 1
%         initRefR = estimateR(refSeq);
%       end
%     end
    
    %% doing innerblock learning to approach the frequency
    [caseSeqAfterInnerBlockLearning] = gInnerSeqLearning(caseSeq, refSeq, blocks, false);
    
    
    %% doing interblock learning starting from new sequence
    [result] = genotypeInterBlockRecomb(caseSeq, caseSeqAfterInnerBlockLearning, blocks, false);
    
    
    %check the sign recover rate
end