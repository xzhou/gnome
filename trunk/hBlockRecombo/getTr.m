function Tr = getTr(Y, r_S, r_R, maskMatrix, signMatrix)
if nargin == 3
    A2 = (2*Y'-1)*(2*Y-1);
    Tr = sum(sum((r_S - r_R).* A2))/2;
else
    r_S = abs(r_S).*signMatrix; %replace the sign
    r_S = r_S.*maskMatrix;
    r_R = r_R.*maskMatrix;
    A2 = (2*Y'-1)*(2*Y-1);
    Tr = sum(sum((r_S - r_R).* A2))/2;
end
end