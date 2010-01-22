function [] = change_env()
%CHANGE_ENV will handle all environmental problems
    path = pwd;
    if start_with(path, '/home/xzhou')
        disp 'running on xzhou.info'
        addpath '~/research_linux/gnome/bioWorkspace/genomeprj/genotypeRecomb'
        addpath '~/research_linux/gnome/bioWorkspace/genomeprj/hBlockRecombo'
        addpath '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis'
        addpath '~/research_linux/gnome/bioWorkspace/genomeprj/readFiles'
        cd '~/research_linux/gnome/bioWorkspace/genomeprj/data/1500DataAnalysis/WTCCC1/TPED'      
    elseif start_with(path, '/u')
        disp 'running on hulk.cs.indiana.edu'
    elseif start_with(path, 'D:')
        disp 'running on Windows, please update th path here :)'
    end
end

function [result] = start_with(str1, str2)
    [m len] = size(str2);
    x = ((str1(1:len) == str2) == 0);
    if sum(x) > 0
        result = 0;
    else
        result = 1;
    end
end