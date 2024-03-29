function [caseSeq4 refSeq4] = randomSelect(sequenceAll, number, refPool)
%Function random divide the sequenceAll into 2 groups and convert it into
%caseSeq and refSeq.
%Argument Number defines the number of people in case group.
%Argument refPool defines how to construct the reference group. 
%0:constuct reference that has the same number with case group.
%1:Use all the remaining people to construct the reference group.


seqAll = sequenceAll;

totalLength = length(seqAll);

if(mod(totalLength, 4)==2)
seqAll(end) = [];
seqAll(end) = [];
end

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));   %reset the random number generator

totalLength = length(seqAll);

if nargin == 1
    numbers = totalLength/4;
else
    numbers = number;
end

sel = [1:totalLength/2];
a = totalLength/2;
b = randi(a,1,1);
selS = zeros(1, totalLength/4);
selR = zeros(1, totalLength/4);

if refPool==0
    for i = 1:totalLength/4
        c = randi(a,1,1);
        selS(1,i) = sel(1, c);
        sel(c) = [];
        a = a-1;
    end
    selR = sel;
    for i=1:numbers
    seqS(i*2-1) = seqAll(selS(i)*2-1);
    seqS(i*2) = seqAll(selS(i)*2);
    
    seqR(i*2-1) = seqAll(selR(i)*2-1);
    seqR(i*2) = seqAll(selR(i)*2);
    end
else
    for i = 1:numbers;
        c = randi(a,1,1);
        selS(1,i) = sel(1, c);
        sel(c) = [];
        a = a-1;
    end
    selR = sel;
    for i=1:numbers
    seqS(i*2-1) = seqAll(selS(i)*2-1);
    seqS(i*2) = seqAll(selS(i)*2);
    end
    for i=1:length(selR)
    seqR(i*2-1) = seqAll(selR(i)*2-1);
    seqR(i*2) = seqAll(selR(i)*2);
    end
end

    

% %for i = 1:totalLength/4
% for i=1:numbers
%     seqS(i*2-1) = seqAll(selS(i)*2-1);
%     seqS(i*2) = seqAll(selS(i)*2);
%     
%     seqR(i*2-1) = seqAll(selR(i)*2-1);
%     seqR(i*2) = seqAll(selR(i)*2);
% end

seqS = seqS';
seqR = seqR';


nS = length(seqS);
Len = length(seqS(1).Sequence);
nR = length(seqR);

int4S = zeros(nS, Len);
int4R = zeros(nR, Len);

for i = 1:nS
    int4S(i,:) = nt2int(seqS(i).Sequence) - 1;
end

for i = 1:nR
    int4R(i,:) = nt2int(seqR(i).Sequence) - 1;
end

caseSeq4 = int4S;
refSeq4 = int4R;

end

