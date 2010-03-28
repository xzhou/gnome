function [caseSeq4 refSeq4 testSeq4] = randomSelectTest(sequenceAll, number)
%Function randomly divides the sequenceAll into 3 groups and convert it into
%caseSeq, refSeq and testSeq.
%Argument Number defines the number of people in case group.
%Argument refPool defines how to construct the reference group. 
%0:constuct reference that has the same number with case group.
%1:Use all the remaining people to construct the reference group.


[nS ~] = size(sequenceAll);

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));   %reset the random number generator

sel = [1:nS];
a = nS;
selS = zeros(1, number);
selR = zeros(1, number);
selT = zeros(1, number);

for i = 1:number
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

for i=1:number
    refSeq4(i,:) = sequenceAll(selR(i), :);
    
    caseSeq4(i,:) = sequenceAll(selS(i), :);
    
    testSeq4(i,:) = sequenceAll(selT(i), :);

end


end

