function scores=computeDrawScores(haps,pop,ratio)

% haps: haplotypes
% pop: population
% ratio: the weight for haplotypes appeared in pop
% scores: scores assigned to each haplotype

if nargin == 2, ratio = 1; end % ratio in default is 1

% assign scores for haplotypes appeared in the population
% the strategy is simple: each haplotype appeared in the population
% receives the same score
[row,col]=size(pop);
scores(length(haps))=0;
for i=1:row
    for j=1:length(haps)
        equal = 1;
        for k=1:col
            if pop(i,k)~=haps(j,k)
                equal = 0;
                break;
            end
        end
        if equal == 1
            scores(j) = 1;
            break;
        end
    end
end

% assign scores for other haplotypes
value = (sum(scores)*(1-ratio)/ratio)/sum(scores==0);
for i=1:length(scores)
    if scores(i)==0
        scores(i) = value;
    end
end
