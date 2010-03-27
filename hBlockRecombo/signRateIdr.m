function [idr, T] = signRateIdr(caseSeq, testSeq, fdr, level, caseR, refR)
%in this experimtns, we try to find the relationship between identification
%rate and sign rate, 

%% get threash hold value for different p level
p = zeros(level, 1);
for i = 1:level
    p(i) = i*level/100;
end
T = findCut(caseR, p);
idr = zeros(1, level);

%% calculate identification rate
for i = 1:level
    Ti = T(i);
    mask = ones(size(caseR));
    mask(abs(caseR)<Ti) = -1;
    caseMaskR = caseR.*mask;
    idr(i) = getIdr(caseSeq, testSeq, fdr, caseMaskR, refR);
end

end


function [T] = findCut(caseR, p)
%given a percentage or multiple percentages, find the value cut of a matrix such that caseR < T has
%p elements
v = upperEle(abs(caseR));
sortv = sort(v, 'descend');
n = length(v);
T = p;
pLen = length(p);
for i = 1:pLen
    pi = p(i);
    idx = round(pi*n);
    T(i) = sortv(idx);
end
end

function [v] = upperEle(M)
%%put the upper triangle elements of M in to v, M must a square matrix
[m ~] = size(M);
index = triu(ones(m), 1);
v = M(logical(index));
end