function [h] = plotPowerDist(m, signMatrix)
    [nrow ncol] = size(m);
    h = figure();
    hold on;
    for i = 1:nrow
        for k = 1:ncol
            x_init = [k-1, k];
            y_init = [i-1, i];
            z_init = double([m(i,k) m(i,k); m(i,k) m(i,k)]);
            box = surface(x_init, y_init, z_init);
        end
    end
    
    maxm = max(max(abs(m)));
    
    load('depColor', 'depColor');
    set(gcf, 'ColorMap', depColor);
    axis equal;
    caxis([-maxm maxm]);
    colorbar;
    
    if nargin == 2
        signVector = MatrixToVec(signMatrix);
        nrow = length(signVector);
        z = ones(nrow, 1) * maxm + 0.1;

        scatter3(signVector(:,1), signVector(:,2), z,'marker', 'x');
    end
    
    h = gcf;
end
