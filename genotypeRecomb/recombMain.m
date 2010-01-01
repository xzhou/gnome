%this is the main function
function [] = recombMain()
    %% add path
    try
        addpath '/home/xzhou/research_linux/gnome/workspace/hBlockRecombo'
        addpath '/home/xzhou/research_linux/gnome/workspace/MyDrTest'
        addpath '/home/xzhou/research_linux/gnome/workspace/genotypeRecomb'
    catch exception
        %do nothing
    end
    
    %% cd to data
    try
        cd 'D:\IUBResearch\Projects\Bioinfor\data\HAPMAP';
        disp 'WINDOWS'
    catch e
        cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP';
        disp 'LINUX'
    end
    
    try
        startParallel(2);
    catch e
        disp "failed to run in parallele"
    end
    
    %% reading genotype data
    fastaFile = 'hapmap_chr7_80SNP_CEU_haplotype.fasta';
    blocks = [1 15; 16 59; 60 77];      %static block assignment
    
    all = readSeq4(fastaFile);
    realR = calcR(all);
    
    genotypeAll = genotypeHelpFuncs.readGenotypeFromFasta(fastaFile);
    
    %% test
    [totalR pA counts] = estimateR(genotypeAll);
%    return;
    
    %% for experiments
    [caseSeq, refSeq] = randomSelectGenotype(genotypeAll);
    
    %calculate target estimated R
    parfor pi = 0:1
      if pi == 0
        targetR = estimateR(caseSeq);
      elseif pi == 1
        initRefR = estimateR(refSeq);
      end
    end
    
    %% doing innerblock learning
    [seq] = gInnerSeqLearning(caseSeq, refSeq, blocks, false);
    
    %% doing interblock learning
    
    
end