function [tp] = getTp(Y, pM1, pM0)
%calculate Tp value using Homer's attack
p = (-(Y-1)+1)/2;
tp = sum(abs(p - pM0) - abs(p - pM1));
end