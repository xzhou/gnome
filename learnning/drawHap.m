function k=drawHap(scores)

% draw a haplotype according to the scores
% scores: each number in scores represents a score of a haplotype

k = 0;
s = sum(scores);
value = s*rand(1);
begin = 0;
for i=1:length(scores)
    e = begin + scores(i);
    if value > begin && value <= e
        k = i;
        break;
    end
    begin = e;
end