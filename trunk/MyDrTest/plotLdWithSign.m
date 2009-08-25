function [] = plotLdWithSign(r, signMatrix, T)
% plotLdWithSign will plot the ld map of the r value and marker the
% different sign
% T is the mask
    
    [numrow, numcol] = size(r);
    
    %filter
    maskIndex = logical(abs(r)<T);
    r(maskIndex) = 0;
    
    figure;
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
    
    load('MyColorMap', 'mycmap');
    set(gcf, 'ColorMap', mycmap);
    axis equal;
    %plot sign
    signVector = MatrixToVec(signMatrix, maskIndex');
    nrow = length(signVector);
    z = ones(nrow, 1);
    %scatter3(signVector(:,1), signVector(:,2), z,'marker', 'x');
end

function [] = plotWithSign(r, T)
% plotLdWithSign will plot the ld map of the r value and marker the
% different sign
% T is the mask
    
    [numrow, numcol] = size(r);
    
    %filter
    maskIndex = logical(abs(r)<T);
    r(maskIndex) = 0;
    
    figure;
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
    
    load('MyColorMap', 'mycmap');
    set(gcf, 'ColorMap', mycmap);
    axis equal;
    %plot sign
    
    signMatrix = sign(r);
    
    signVector = MatrixToVec(signMatrix, maskIndex');
    nrow = length(signVector);
    z = ones(nrow, 1);
    %scatter3(signVector(:,1), signVector(:,2), z,'marker', 'x');
end