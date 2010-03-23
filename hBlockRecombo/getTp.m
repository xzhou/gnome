function [tp] = getTp(Y, pM0, pM1)
%calculate Tp value using Homer's attack
tp = sum(abs(Y - pM0) - abs(Y - pM1));
end