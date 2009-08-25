function [x] = MatrixToVec(signMatrix, maskIndex)
    if nargin == 0;
        load;
        signMatrix = signDiff;
    elseif nargin == 1
        [m,n] = size(signMatrix);
        maskIndex = zeros(m,n);
    end
    

    [m n] = size(signMatrix);
    
    x = zeros(0,2);
    for i = 1:m
        for j = 1:n
            if(signMatrix(i,j) == -1 && maskIndex(i,j) ~= 1)
                x = [x; i-0.5 j-0.5];
            end
        end
    end
    
    %scatter(x(:,1), x(:,2), 'marker', 'x');
    
end