function [] = innerBlockLearning()

cd 'D:\IUBResearch\Projects\Bioinfor\data\HAPMAP';
 
    %% reading genotype data
    fastaFile = 'hapmap_chr7_80SNP_CEU_haplotype.fasta';
    blocks = [1 15; 16 59; 60 77];      %static block assignment
    
    all = readSeq4(fastaFile);
    realR = calcR(all);
    
    genotypeAll = genotypeHelpFuncs.readGenotypeFromFasta(fastaFile);
    
    %% init r
    [totalR pA counts] = estimateR(genotypeAll);
    %%return;
    
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
    
    type innerBlockMutate.m
    
    type innerBlockFitness.m
    
    
    fitnessfcn = @(x) innerBlockFitness(x,caseSeq(:,1:15));
    
    options = saoptimset('DataType', 'custom', 'AnnealingFcn', @innerBlockMutate, ...
    'StallIterLimit',800, 'ReannealInterval', 800, 'PlotInterval', 50);
    %%
    % Finally, we call simulated annealing with our problem information.
    genotypeBlock = simulannealbnd(fitnessfcn, refSeq(:,1:15), [], [], options);
    
    test =  magic(6);
    
    