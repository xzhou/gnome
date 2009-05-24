function [s,err]=computeSingleSNPFrequencyError(A,Aref,first,last)

% compute the error of single SNP frequencies between two populations A and
% Aref.
% A: individual sequences
% Aref: reference sequences
% first: index of the first SNP
% last: index of the last SNP
% s: sum of the errors
% err: error

[row,col]=size(A);
if nargin <= 2
    first = 1;
    last = col;
end

err(last-first+1) = 0;
for i=first:last
    err(i) = abs(sum(A(:,i))-sum(Aref(:,i)))/row;
end
s = sum(err)/length(err);