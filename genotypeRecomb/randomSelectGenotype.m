function [caseGenotype, refGenotype, caseID, refID] = randomSelectGenotype(genotype, idInfo, number)
%%RANDOM_SELECT_GENOTYPE will randomly select two groups of people from genotype    
    seqAll = genotype;
    [rows columns] = size(seqAll);
    
    if 2*number > rows
        MException('GenotypeRecomb:randomSelectGenotype', 'not enough individuals');
    end

    if(mod(rows, 2) == 1)
        seqAll(end, :) = [];
    end

    RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

    [nIndividuals nSnps] = size(seqAll);

    %randomly select the case and reference group
    rowIndexs = 1:nIndividuals;
    a = nIndividuals;
    selectedID = [];
    for i = 1:(2*number)
       r = randi(a, 1, 1);
       selectedID = [selectedID, rowIndexs(r)];
       rowIndexs(r) = [];
       a = a - 1;
    end
    
    %randomly select case from the randomly selected population
    caseIndex = [];
    a = 2*number;
    for i = 1:number
       c = randi(a, 1, 1);
       caseIndex = [caseIndex, selectedID(c)];
       selectedID(c) = [];
       a = a - 1;
    end
    
    referenceIndex = selectedID;
    
    caseGenotype = genotype(caseIndex, :);
    refGenotype = genotype(referenceIndex, :);
    
    caseID = idInfo(caseIndex);
    refID  = idInfo(referenceIndex);
end