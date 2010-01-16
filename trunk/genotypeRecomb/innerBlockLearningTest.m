function [] = innerBlockLearning()
%This is a function for testing the SA algorithm

cd 'D:\IUBResearch\Projects\Bioinfor\data\HAPMAP';

%% reading genotype data
fastaFile = 'hapmap_chr7_80SNP_CEU_haplotype.fasta';
blocks = [1 15; 16 59; 60 77];      %static block assignment
[m n] = size(blocks);

all = readSeq4(fastaFile);
realR = calcR(all);

genotypeAll = genotypeHelpFuncs.readGenotypeFromFasta(fastaFile);

%% init r
[totalR pA counts] = estimateR(genotypeAll);
%%return;

%% for experiments
%[caseSeq, refSeq] = randomSelectGenotype(genotypeAll);

[caseSeq, refSeq] = randomSelect(all);

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

seqAfterInnerLearning = zeros(size(caseSeq));

%Use For for each block
for i = 1:m
    %fitnessfcn = @(x) innerBlockFitness(x,caseSeq(:,1:15));
    fitnessfcn = @(x) innerBlockFitness(x,caseSeq(:,blocks(i,1):blocks(i,2)));
    
    options = saoptimset('DataType', 'custom', 'AnnealingFcn', @innerBlockMutate, ...
        'StallIterLimit',800, 'ReannealInterval', 800);
    % Finally, we call simulated annealing with our problem information.
    genotypeBlock = simulannealbnd(fitnessfcn, refSeq(:,blocks(i,1):blocks(i,2)), [], [], options);
    %genotypeBlock = simulannealbnd(fitnessfcn, refSeq, [], [], options);
    
    seqAfterInnerLearning(:,blocks(i,1):blocks(i,2)) = genotypeBlock;
end


%%May add the start point caculation process here.

%%Start Inter Block Learning.

type interBlockMutate.m

type interBlockFitness.m

signMatrix = zeros(size(caseSeq));
signMatrix = sign(estimateR(seqAfterInnerLearning));


for i = 1:m
    blocks(i, 3) = blocks(i,2)-blocks(i,1)+1;
end
for i= 1:(m-1)
    for j = i+1:m
        block1 = blocks(i,:);
        block2 = blocks(j,:);
        if(block1(1,3)>= block2(1,3))
            fitnessfcn = @(x) interBlockFitness(x,caseSeq,block1,block2);
            %interBlockMutateNew = @(optimValues, problemData) interBlockMutate(optimValues, problemData, block2);
            options = saoptimset('DataType', 'custom', 'AnnealingFcn', @interBlockMutate, ...
                'StallIterLimit',800, 'ReannealInterval', 800);
            genotypeBlock = simulannealbnd(fitnessfcn, seqAfterInnerLearning, [], [], options);
            
        else
            fitnessfcn = @(x) interBlockFitness(x,caseSeq,block1,block2);
            %interBlockMutateNew = @(optimValues, problemData) interBlockMutate(optimValues, problemData, block1);
            options = saoptimset('DataType', 'custom', 'AnnealingFcn', @interBlockMutate, ...
                'StallIterLimit',800, 'ReannealInterval', 800);
            genotypeBlock = simulannealbnd(fitnessfcn, seqAfterInnerLearning, [], [], options);
        end
    end
    %Get the sign from the seq after learning at the end
    blockMask = getBlockMaskForEval(caseSeq, block1, block2);
    signMatrix = signMatrix.*(blockMask==0) + sign(estimateR(genotypeBlock)).*blockMask;
    
end

%Combine the sign results get from inner and inter learning.

test =  magic(6);

end

