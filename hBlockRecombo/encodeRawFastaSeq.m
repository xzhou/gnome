function [hap01Seq, alleleMapping] = encodeRawFastaSeq(rawFastaSeq, alleleMapping)
%convert raw fasta sequence to 01 encoding, 1 is major, 0 is minor
hapIntSeq = seq2int(rawFastaSeq);
hapSeqNoID = getSeqMatrix(hapIntSeq);
if nargin == 1
    alleleMapping = getMajorAllele(hapSeqNoID);
end

[m n] = size(hapSeqNoID);
hap01Seq = zeros(m, n);
for i = 1:m
    hap01Seq(i, :) = (hapSeqNoID(i,:) == alleleMapping) + 0;
end
end