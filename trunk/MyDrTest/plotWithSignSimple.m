function [h] = plotWithSignSimple(r1, r2, T)
% r1 is sample, r2 is case
% T is the mask threshold, 
    
    if nargin == 1
        r2 = r1;
        T = 0.0;
    elseif nargin == 2
        T = 0.0;
    end
    
    [numrow, numcol] = size(r1);
    
    %mix r value
    r = r1;
    index = ones(numrow);
    index = logical(triu(index));
    r = r2;
    r(index) = r1(index);
    
    %filter
    maskIndex = logical(abs(r)<T);
    r(maskIndex) = 0;
    
    h = figure();
    zdir = [0 0 1];
    center = [numrow/2 numrow/2 0];
    rotate(h,zdir,45,center);
    hold on;
    %plotLD map
    for i = 1:numrow
        for k = 1:numcol
            x_init = [k-1, k];
            y_init = [i-1, i];
            z_init = double([r(i,k) r(i,k); r(i,k) r(i,k)]);
            box = surface(x_init, y_init, z_init);
        end
    end
    
    load('depColor', 'depColor');
    set(gcf, 'ColorMap', depColor);
    axis equal;
    %plot sign
    
    signMatrix = sign(r1).*sign(r2);
    
    signVector = MatrixToVec(signMatrix, maskIndex');
    nrow = length(signVector);
    z = ones(nrow, 1);
    scatter3(signVector(:,1), signVector(:,2), z,'marker', 'x');
    colorbar;
    caxis([-1 1]);
    h = gcf;
end