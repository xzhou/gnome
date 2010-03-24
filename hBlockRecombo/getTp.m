function [tp] = getTp(Y, pM0, pM1)
%calculate Tp value using Homer's attack
p = (-(Y-1)+1)/2;
tp = sum(abs(p - pM1) - abs(p - pM0));
end