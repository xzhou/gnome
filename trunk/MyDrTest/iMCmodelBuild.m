function [iMCmodel sim] = iMCmodelBuild(fastafile, pseudocount)
% building inhomogenous Markov chain from sample;
% sample is an m*n marix with m=(number of sequences), n=(length of each
% sequence)
% 
if nargin < 2, pseudocount = 0; end
if nargin < 1, fastafile = 'chr10_FGFR2_200kb_phased.fasta'; end;
if strcmp(class(fastafile),'char') | strcmp(class(fastafile),'struct')
    % the input is the name of the fasta file, or a struct with sequences
    if strcmp(class(fastafile),'char')
        seq = fastaread(fastafile);
    else
        seq = fastafile;
    end
    nSample = length(seq);
    Len = length(seq(1).Sequence);
    nState = 4;

    sample = zeros(nSample,Len);
    for i = 1:length(seq)
        sample(i,:)=nt2int(seq(i).Sequence) - 1;
    end
elseif strcmp(class(fastafile),'double')
    % the input sequences have been converted to int already
    sample = fastafile;
    [nSample, Len] = size(sample)
    nState = 4;    
end

sim = zeros(nSample,nSample);

for i = 1:nSample
    for j = i+1:nSample
        sim(i,j) = sum(sample(i,1:Len) == sample(j,1:Len));
    end
end
sim = sim;

% tmp = tabulate(sample(:,1));
for s = 1:nState
    iMCmodel.initial(s) = sum(sample(:,1)==s-1);
end
iMCmodel.initial = iMCmodel.initial + pseudocount;
iMCmodel.initial = iMCmodel.initial/sum(iMCmodel.initial);

iMCmodel.transition = zeros(nState, nState, Len-1);
for i = 2:Len
    % [tmp c p l] = crosstab(sample(:,i-1),sample(:,i)) + pseudocount;
    for s1 = 1:nState
        for s2 = 1:nState
            tmp(s1,s2)=sum(sample(:,i-1)==s1-1 & sample(:,i)==s2-1);
        end
    end
    tmp = tmp + pseudocount;
    iMCmodel.transition(:,:,i-1) = tmp./repmat(sum(tmp,2),1,size(tmp,2));
end
