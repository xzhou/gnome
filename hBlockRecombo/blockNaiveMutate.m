function [currentFreq] = blockNaiveMutate(currentFreq)
%naive mutate, each time randomly select two hyplotyp frequency and change
%them
    [m, n] = size(currentFreq);
    if n ~= 1
        e = MException('blockNaiveMutate:column', 'frequency must by n by 1');
        throw(e);
    end
    x1 = randi(m);
    x2 = randi(m);
    while x1 == x2 || currentFreq(x2, 1) == 0
        x1 = randi(m);
        x2 = randi(m);
    end
    
    currentFreq(x1, 1) = currentFreq(x1, 1) + 1;
    currentFreq(x2, 1) = currentFreq(x2, 1) -1;
end