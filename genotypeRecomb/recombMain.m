%this is the main function
function [] = recombMain()
    
    

    startParallel(2);
    
    %% reading genotype data
    fastaFile = 'hapmap_chr7_80SNP_CEU_haplotype.fasta';
    blocks = [1 15; 16 59; 60 77];      %static block assignment
    
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