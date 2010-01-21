%this function divide the genotype sequence into two groups
function [caseGenotype, refGenotype] = randomSelectGenotype(genotype)

seqAll = genotype;
[rows columns] = size(seqAll);

if(mod(rows, 2) == 1)
    seqAll(end, :) = [];
end

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

[rows columns] = size(seqAll);

sel = [1:rows];
a = rows;
selS = zeros(1, rows/2);
selR = zeros(1, rows/2);

for i=1:rows/2
    c=randi(a,1,1);
    selS(1,i) = sel(1,c);
    sel(c) = [];
    a = a-1;
end

selR = sel;

for i=1:rows/2
    caseGenotype(i,:) = seqAll(selS(i),:);
    refGenotype(i,:) = seqAll(selR(i),:);
end

end