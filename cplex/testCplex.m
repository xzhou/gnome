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
    
    rhs = 9.3;
    
    cplex.Model.name = 'test mip';
    cplex.Model.lb = lb;
    cplex.Model.ub = ub;
    cplex.Model.A = A;
    cplex.Model.rhs = rhs;
    
    cplex.Model.ctype = 'II';
    cplex.Model.sense = 'maximize';
    cplex.Param.mip.pool.relgap.Cur = 0.1; 
    cplex.populate();
    
end
