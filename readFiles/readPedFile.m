function [genotype, alleleMapping, idInfo] = readPedFile(fileName)
%READ_PED_FILE will read genotype sequence from a ped file and save it to a
%matrix
    if nargin == 0
        M = importdata('Affx_gt_58C_Chiamo_07.tped.extract.inp.ped');
    else
        M = importdata(fileName);
    end
    
    [m n] = size(M);
    
    startParallel(2);
    
    idInfo = [];
    genotype = [];
    parfor i=1:m
        aSeq = M{i, :};
        tokens = strsplit(' ', aSeq);
        idInfo = [idInfo, tokens(:, 2)]
        tokens = tokens(:, 7:end);
        tokens = [tokens{1:end}];
        genotype = [genotype; nt2int(tokens)];
    end
    
    %get major allele
    [m n] = size(genotype);
    majorAllele = zeros(1, n/2);
    for j = 1:(n/2)
        seqA = genotype(:, 2*j -1);
        seqB = genotype(:, 2*j);
        seq = [seqA', seqB'];
        cA = sum(seq==1);
        cC = sum(seq==2);
        cG = sum(seq==3);
        cT = sum(seq==4);
        counts = [cA, cC, cG, cT];
        [a maxIdx] = max(counts);
        majorAllele(1, j) = maxIdx;
    end
    
    newGenotype = zeros(m, n/2);
    for i = 1:m
        for j = 1:n/2
            A = genotype(i, 2*j-1);
            B = genotype(i, 2*j);
            C = majorAllele(1, j);
            if A == C && B == C
                newGenotype(i,j) = 0;
            elseif A~= C && B~=C
                newGenotype(i,j) = 2;
            else
                newGenotype(i,j) = 1;
            end
        end
    end
    
    %return value
    genotype = newGenotype;
    alleleMapping = majorAllele;
    
end