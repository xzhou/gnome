function [caseSeq4 refSeq4 testSeq4] = randomSelectTestNew(sequenceAll, number)
%Function random divide the sequenceAll into 2 groups and convert it into
%caseSeq and refSeq.
%Argument Number defines the number of people in case group.
%Argument refPool defines how to construct the reference group. 
%0:constuct reference that has the same number with case group.
%1:Use all the remaining people to construct the reference group.


seqAll = sequenceAll;

totalLength = length(seqAll);

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));   %reset the random number generator

totalLength = length(seqAll);

if nargin == 1
    numbers = totalLength/6;
else
    numbers = number;
end

sel = [1:totalLength/2];
a = totalLength/2;
selS = zeros(1, numbers);
selR = zeros(1, numbers);
selT = zeros(1, numbers);

for i = 1:numbers
    r = randi(a,1,1);
    selR(1,i) = sel(1, r);
    sel(r) = [];
    a = a-1;
    
    c = randi(a,1,1);
    selS(1,i) = sel(1, c);
    sel(c) = [];
    a = a-1;
    
    t = randi(a,1,1);
    selT(1,i) = sel(1, t);
    sel(t) = [];
    a = a-1;
end

for i=1:numbers
    seqS(i*2-1) = seqAll(selS(i)*2-1);
    seqS(i*2) = seqAll(selS(i)*2);
    
    seqR(i*2-1) = seqAll(selR(i)*2-1);
    seqR(i*2) = seqAll(selR(i)*2);
    
    seqT(i*2-1) = seqAll(selT(i)*2-1);
    seqT(i*2) = seqAll(selT(i)*2);
end


seqS = seqS';
seqR = seqR';
seqT = seqT';


nS = length(seqS);
Len = length(seqS(1).Sequence);

int4S = zeros(nS, Len);
int4R = zeros(nS, Len);
int4T = zeros(nS, Len);

for i = 1:nS
    int4S(i,:) = nt2int(seqS(i).Sequence) - 1;
    
    int4R(i,:) = nt2int(seqR(i).Sequence) - 1;
    
    int4T(i,:) = nt2int(seqT(i).Sequence) - 1;
end


caseSeq4 = int4S;
refSeq4 = int4R;
testSeq4 = int4T;

end

