function testCplex()

%testing how to get multiple solutions
% target: max x + y
% subject to:
%   x <= 5.5
%   y <= 5.5
%   x + y <= 9.3
% 

    cplex = Cplex();
    
    obj = [1, 1]';
    A = [1, 1];
    
    lb = [0, 0]';
    ub = [5.5, 5.5]';
    
    rhs = [9.3];
    
    cplex.Model.name = 'test mip';
    cplex.Model.lb = lb;
    cplex.Model.obj = obj;
    cplex.Model.ub = ub;
    cplex.Model.A = A;
    cplex.Model.rhs = rhs;
    cplex.Model.lhs = [-inf];
    cplex.Model.sense = 'maximize';
    
    %optional
    cplex.Model.ctype = 'II';
    cplex.Param.mip.pool.intensity.Cur = 4;
    %cplex.Param.mip.limits.populate.Cur = 3;
    cplex.Param.mip.pool.relgap.Cur = 0;
    cplex.populate();
    %control searching strategy
    numsoln = size (cplex.Solution.pool.solution, 1);
    
    for i = 1:numsoln
      %fprintf ('%8d  %13g', i, cplex.Solution.pool.solution(i).objval);
      cplex.Solution.pool.solution(i).x
   end
end
