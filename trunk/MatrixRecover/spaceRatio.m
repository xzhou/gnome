function [r] = spaceRatio(N, L)
%calculate the space ratio
%we use nchoose k to represent the space for large L, but for samll L, we
%need exact space
S = nchoosek(2^L, N);
D = (N+1)^(nchoosek(L,2)+1);
r = S/D;
end