function genotype2Blocks = interBlockMutate(optimValues, problemData, block1)
%function genotype2Blocks = interBlockMutate(optimValues, problemData,
%caseBlock, block1, block2)
%This function defines the interBlock mutate algorithm.

genotype2Blocks = optimValues.x;
%The whole current ref block seq.
for i=1:floor(optimValues.temperature)+1
[m n] = size(genotype2Blocks);
genotype2Blocks = neighbor(genotype2Blocks,m,n);
end

function genotype2Blocks = neighbor(genotype2Block, m, n)
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));   %reset the random number generator

idx1 = randi([1 m]);
idx2 = randi([1 m]);

a = block1(1, 1);
b = block1(1, 2);

while genotype2Blocks(idx1, a:b) == genotype2Blocks(idx2, a:b)
    idx1 = randi([1 m]);
    idx2 = randi([1 m]);
end

%swap the sequence
temp = genotype2Blocks(idx1, a:b);
genotype2Blocks(idx1, a:b) = genotype2Blocks(idx2, a:b);
genotype2Blocks(idx2, a:b) = temp;
