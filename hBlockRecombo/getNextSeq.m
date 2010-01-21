function [newSeq, newBlockFreq] = getNextSeq(hyplotypes, currentBlockFreq)
    newBlockFreq = blockNaiveMutate(currentBlockFreq);
    newSeq = blockReconstruct(hyplotypes, newBlockFreq);
end