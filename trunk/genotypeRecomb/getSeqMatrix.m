function [intSeqMatrix] = getSeqMatrix(intSeqWithID)
%remove the id information from sequence to and make it a matrix for easy
%processing
    [tmp, nSeq] = size(intSeqWithID);
    [seqLen] = length(intSeqWithID(1).Sequence);
    intSeqMatrix = zeros(nSeq, seqLen);
    for i = 1:nSeq
        intSeqMatrix(i,:) = intSeqWithID(i).Sequence;
    end
end