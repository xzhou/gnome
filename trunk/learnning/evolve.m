function evolve(initPop,finalPop,realPop,hapPop,blksize,fun,first,last,firstSNP,lastSNP)

% evolve('initHap100.mat','finalHap100.mat','cse100.mat','haps800.mat',5,@dp
% rime,1,100,1,80)
% Function
% Starts from the initial haplotype population, evolve to minimize the error
% Arguments
% initPop: initial population initHap100.mat
% finalPop: final population, next population finalHap100.mat
% realPop: real population cse100.mat
% hapPop: Hapmap population haps800.mat
% blksize: the size of each block 5
% first: the index of the first individual 
% last: the index of the last individual
% firstSNP: the index of the first SNP
% lastSNP: the index of the last SNP
%load(realPop,'cse'); % load variable cse
%load(hapPop,'ref');  % load variable ref
%load(initPop, 'evolver');       % load variable evolver

cse = importdata(realPop);
ref = importdata(hapPop);
evolver = importdata(initPop);

cse = cse(first:last,firstSNP:lastSNP); % real case population
evolver = evolver(first:last,firstSNP:lastSNP); %current population
ref = ref(:,firstSNP:lastSNP);  %hapref population

rs = computeRSquaresForBlock(cse,fun); % compute rsquares from real population

[row,col]=size(cse);

haps = generateHaplotypes(blksize);
ratio = 0.7;  % The ratio for haplotypes appeared in Hapmap population. For details, please refer to computeDrawRate(...)
first = firstSNP;
col = lastSNP;
index = 1;
while first <= col
    last = first+blksize-1;
    if last > col, last = col; end

    scores(index,:)=computeDrawScores(haps,ref(:,first:last),ratio);
    blks(index,:)=[first last]; %block is deliminator position

    first = last+1;
    index = index + 1;
end

tol = 1e-14; % tolerance
scale = 0.000001; % scale
num = 1; % # of replaced records for each generation
seeds = 5; % the # of trials for each generation

r = computeRSquaresForBlock(evolver,fun); % compute rsquares from evolving population

rerr = computePairwiseFrequencyError(rs,r,firstSNP,lastSNP)*scale;  %parewise error
perr = computeSingleSNPFrequencyError(evolver,cse,firstSNP,lastSNP); %single snp frequency error
L = (perr+rerr); % error

select(seeds) = 0;
err(seeds) = 0;
while L > tol % continue if the error is larger than tol
    orig = evolver; 
    [blkRow,blkCol]=size(blks);
    for i = 1:blkRow
        first = blks(i,1);
        last  = blks(i,2);
        b = firstSNP;
        e = lastSNP;
        if b < firstSNP, b=firstSNP; end
        if e > lastSNP, e=lastSNP; end
        for k=1:num
            iHap=drawHap(scores(i,:)); % choose a haplotype
            for n=1:seeds
                select(n)=randInt(1,row); % randomly select a record
                if n ~= 1
                    x = select(n); %randomly select a row
                    if last - first == blksize - 1
                        evolver(x,first:last)=haps(iHap,:); % replace the selected record with the chosen haplotype
                    else
                        tmpHap = haps(iHap,:);
                        evolver(x,first:last) = tmpHap(1:last-first+1);
                    end
                    r = updateBlock(evolver,r,fun,first,last); % evolver is changed, so the corresponding r should be updated
                end
                rerr = computePairwiseFrequencyError(rs,r,b,e)*scale;
                perr = computeSingleSNPFrequencyError(evolver,cse,b,e);
                err(n) = perr+rerr;
                if n ~= 1, 
                    evolver(select(n),first:last)=orig(select(n),first:last);
                end
            end
            [value,index]=min(err);
            if index ~=1
                randomRow = select(index);
                if last - first == blksize - 1
                    evolver(randomRow,first:last)=haps(iHap,:); % update evolver
                else
                    tmp = haps(iHap,:);
                    evolver(randomRow, first:last) = tmp(1:last-first+1);
                end
                save(finalPop,'evolver');  % save to file
                rate = recoverRate_ex(evolver,cse); % recovery rate
                fprintf(1,'blk %d improved,value=%.6f,orig=%.6f,recovered=%.4f\n',i,value,err(1),rate);
            end
            r = updateBlock(evolver,r,fun,first,last);
        end
    end
end