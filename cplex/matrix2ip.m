function matrix2ip(M)
if nargin == 0
    M = [0 0 0;0 1 1; 1 0 1; 1 1 0; 1 1 0];
    M = [0 0 0 1 1;0 1 1 0 0; 1 0 1 0 1; 1 1 0 0 0; 1 1 0 0 1]
end
 
%M = M(1:2, 1:end);
    
[m n] = size(M);

%single allele frequency
p = sum(M)';

%calculate pairwise allele frequency
c00 = matrix2c00(M);

%construct Integer Programming problem
nPair = nchoosek(n,2);  %the number of pairs of pairwise counts
nAuxVar = m*nPair;      %the number of auxiliary variables, each pair of columns has m pair
nVar = m*n + nAuxVar;   %the number of variables including the Aux variables

%target as an feasibility problem
f = zeros(nVar, 1);

Aeq = zeros(n+nPair, nVar);   % n P constraints, 
% p constraints
P = zeros(n, nVar);
x=ones(1,m);
for i = 1:n
    for j = 1:n
        P(i, (m*(i-1)+1):(m*i)) = x;
    end
end

Aeq(1:n, :) = P;
beq = p;

% c00 constraints
CC = zeros(nPair, nAuxVar);
for i = 1:nPair
   b = m*(i-1)+1;
   e = m*i;
   CC(i, b:e) = ones(1, m);
end
Aeq((n+1):end, (m*n+1):end) = CC;
beq = [beq; c00];

lhs = beq;

% bound Aux variables

%AuxA = zeros(2*m*nPairs, nVar);

Aineq = [];
bineq = [];

for j = 1:(n-1)
    for k = (j+1):n
        for i = 1:m             %each pair has m values
            d1 = (j-1)*m+i;     %index of snp j, individual i
            d2 = (k-1)*m+i;     %index of snp k, individual i
            dd_relative = calcdd(i,j,k,m,n);
            dd = dd_relative + m*n;
            tmp = zeros(1, nVar);
            tmp(1, [d1, d2, dd]) = [-1, -1, 2];
            bineq= [bineq; 0];
            Aineq = [Aineq; tmp];   %constraints for Aux variables
            tmp(1, [d1, d2, dd]) = [1, 1, -2];
            Aineq = [Aineq; tmp];   %constraints for Aux variables
            bineq = [bineq; 1]; 
        end
    end
end


A = [Aeq; Aineq];
lhs = [lhs; ones(size(bineq))*-inf];
rhs = [beq; bineq]; %relax constraints but does not change the result

lb = zeros(nVar, 1);    %binary constraints >= 0
ub = ones(nVar, 1);     %binary constrinats <= 1
ctype = [];
for i = 1:nVar
    ctype = [ctype, 'B'];
end

%options = cplexoptimset;
%options.DIagnostics = 'on';
%options.SolnPoolPop = 2;

%formalize cplex solver
%[x, fval, exitflag, output] = cplexmilp ...
    %(f, Aineq, bineq, Aeq, beq, [], [], [], lb, ub, ctype, [], options);

%reversedMatrix = reshape(x(1:(m*n)), m, n)

%build mip model
cplex = Cplex();

cplex.Model.obj = f;
cplex.Model.A = A;
cplex.Model.rhs = rhs;
cplex.Model.lhs = lhs;
cplex.Model.name = 'reverse';
cplex.Model.lb = lb;
cplex.Model.ub = ub;
cplex.Model.sense = 'maximize';
cplex.Model.ctype = ctype;
cplex.writeModel('reverse.lp');

%options
cplex.Param.mip.pool.intensity.Cur = 4;

cplex.populate();

numsoln = size (cplex.Solution.pool.solution, 1)

for i = 1:numsoln
    fprintf(1, 'solutions %d\n', i);
    %cplex.Solution.pool.solution(i).x
    reverseM = reshape(cplex.Solution.pool.solution(i).x(1:(m*n)), m, n)
end

end

%% function calculate the Aux var index
function [x] = calcdd(i,j,k, m, n)
% calculate the location of Aux varp
pairIndex = calcPairIndex(j,k,n);
x = (pairIndex-1) * m + i;
end

