function [ratio] = spaceRatioExact(N, L)
%we calcuate the exact space ratio
matrixSpace = 2^(N*L)/factorial(N);
constraintSpace = (N+1)^(nchoosek(L,2)+1)
ratio = matrixSpace / constraintSpace;
end