function [x] = calcPairIndex( i, j, n)
    if i > j
        x = 0;
        warning('calc:x','error');
        return;
    end
    x = (i-1)*(2*n-i)/2;
    
    x = x + j - i;
end