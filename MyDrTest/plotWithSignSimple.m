function [h] = plotWithSignSimple(r1, r2, T, plotSign)
% r1 is sample, r2 is case
% T is the mask threshold, 
% r1 is the upleft, r2 is lowerright triangle
    
    if nargin == 1
        r2 = r1;
        T = 0.0;
        plotSign = 0;
    elseif nargin == 2
        T = 0.0;
        plotSign = 1;
    elseif nargin == 3
        plotSign = 1;
    end
    
    r1 = double(r1);
    r2 = double(r2);
    
    [numrow, numcol] = size(r1);
    
    %mix r value
    r = r1;
    index = ones(numrow);
    index = logical(triu(index));
    r = r2;
    r(index) = r1(index);
    
    zLim = max(max(abs(r)));
    
    %filter
    maskIndex = logical(abs(r)<T);
    r(maskIndex) = 0;
    
    h = figure();
    zdir = [0 0 1];
    center = [numrow/2 numrow/2 0];
    rotate(h,zdir,45,center);
    xlim([0 numcol]);
    ylim([0 numcol]);
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
    
    %load('depColor', 'depColor');
    %set(gcf, 'ColorMap', depColor);
    axis equal;
    %plot sign
    if plotSign
        signMatrix = sign(r1).*sign(r2);
        signVector = MatrixToVec(signMatrix, maskIndex');
        nrow = length(signVector);
        z = ones(nrow, 1);
        scatter3(signVector(:,1), signVector(:,2), z,'marker', 'x', 'LineWidth', 0.1);
    end
    colorbar;
    caxis([-zLim zLim]);

    h = gcf;
end