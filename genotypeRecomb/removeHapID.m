function [noIDSeq] = removeHapID(haplotypeSeq)
nSeq = length(haplotypeSeq);
noIDSeq = repmat(haplotypeSeq(1).Sequence, nSeq, 1);
for i = 1:nSeq
    noIDSeq(i,:) = haplotypeSeq(i).Sequence;
end
end