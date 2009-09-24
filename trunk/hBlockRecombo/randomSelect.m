function [caseSeq4 refSeq4] = randomSelect(sequenceAll)
%function random divide the sequenceAll into 2 groups and convert it into a
%caseSeq and refSeq

seqAll = sequenceAll;

totalLength = length(seqAll);

if(mod(totalLength, 4)==2)
seqAll(end) = [];
seqAll(end) = [];
end

totalLength = length(seqAll);

sel = [1:totalLength/2];
a = totalLength/2;
b = randi(a,1,1);
selS = zeros(1, totalLength/4);
selR = zeros(1, totalLength/4);

for i = 1:totalLength/4
    c = randi(a,1,1);
    selS(1,i) = sel(1, c);
    sel(c) = [];
    a = a-1;
end

selR = sel;

for i = 1:totalLength/4
    seqS(i*2-1) = seqAll(selS(i)*2-1);
    seqS(i*2) = seqAll(selS(i)*2);
    
    seqR(i*2-1) = seqAll(selR(i)*2-1);
    seqR(i*2) = seqAll(selR(i)*2);
end

seqS = seqS';
seqR = seqR';

nS = length(seqS);
Len = length(seqS(1).Sequence);
nR = length(seqR);

int4S = zeros(nS, Len);
int4R = int4S;

for i = 1:length(seqS)
    int4S(i,:) = nt2int(seqS(i).Sequence) - 1;
    int4R(i,:) = nt2int(seqR(i).Sequence) - 1;
end

caseSeq4 = int4S;
refSeq4 = int4R;

end

