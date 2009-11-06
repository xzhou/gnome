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
    
    genotypeAll = genotypeHelpFuncs.readGenotypeFromFasta(fastaFile);
    
    %% test
    [totalR pA counts] = estimateR(genotypeAll);
%    return;
    
    %% for experiments
    [caseSeq, refSeq] = randomSelectGenotype(genotypeAll);
    
    %%calculate target estimated R
    targetR = estimateR(caseSeq);
    initRefR = estimateR(refSeq);
end