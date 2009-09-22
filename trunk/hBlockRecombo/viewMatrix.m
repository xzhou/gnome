function [g] = viewMatrix(m)
    [numrow, numcol] = size(m);
    h = figure();
    zdir = [0 0 1];
    for i = 1:numrow
        for k = 1:numcol
            x_init = [k-1, k];
            y_init = [i-1, i];
            z_init = double([m(i,k) m(i,k); m(i,k) m(i,k)]);
            box = surface(x_init, y_init, z_init);
        end
    end
    colorbar;
    axis equal; 
    g = gcf;
end