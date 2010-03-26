function [iMCmodel sim] = iMCmodelBuild_N_state(sample, pseudocount)
% building inhomogenous Markov chain from sample;
% sample is an m*n marix with m=(number of sequences), n=(length of each
% sequence)
% 
if nargin < 2, pseudocount = 0; end

nState = length(unique(sample));

if nState ~= 2
    e = MException('error:iMCmodelBuild', 'only 2 type supported');
    throw(e);
end

[nSample, Len] = size(sample);

sim = zeros(nSample,nSample);

% for i = 1:nSample
%     for j = i+1:nSample
%         sim(i,j) = sum(sample(i,1:Len) == sample(j,1:Len));
%     end
% end
% sim = sim;

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
    tmpTrans =  tmp./repmat(sum(tmp,2),1,size(tmp,2));
    tmpTrans(isnan(tmpTrans)) = 0;
    iMCmodel.transition(:,:,i-1) = tmpTrans;
end
