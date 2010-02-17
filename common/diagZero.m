function [x] = diagZero(M)
[m n] = size(M);

x = M;
if m~=n
    warning('M is not consistent');
else
    for i = 1:m
        x(i,i) = 0;
    end
end
end