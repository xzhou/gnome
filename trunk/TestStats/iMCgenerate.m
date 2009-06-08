function sample = iMCgenerate(iMCmodel,sampleSize)
% generate in homogenous morkov chain 
Len = size(iMCmodel.transition,3)+1; % length of sequence
nState = length(iMCmodel.initial); % number of states
sample = zeros(sampleSize,Len);

sample(:,1) = randbyProb(iMCmodel.initial,sampleSize);
for i = 2:Len
    for K = 0:nState-1
        indxK = find(sample(:,i-1) == K);
        numK = length(indxK);
        if numK == 0
            continue
        end
        sample(indxK,i) = randbyProb(iMCmodel.transition(K+1,:,i-1),numK);
    end
end
    
function r = randbyProb(p,n)
% generate r from {0,1, s-1} with probability p = [p1 p2 p3 ... ps];
% n is the number of time to sample (i.e. size(c)==[n,1])
% n = 2;
% p = [0.6 0.1 0.3];
nState = length(p);
for i = 2:nState
    p(i)=p(i)+p(i-1);
end
if p(end)~=1
    error('Not probability vector');
end

% 
r = rand(n,1);
% r
[tmp indx]= sort([p'; r]); % p(j-1)<=r<p(j) will be assigned to class j
IX = find(indx<=nState); % index in indx for the probabilities
IX = [0; IX];
for i = 2:length(IX)
    indx_for_one_class = indx(IX(i-1)+1:IX(i)-1) - nState;
    r(indx_for_one_class)= i-2; % index for random numbers to be assigned to class i-1
end
% r


return;