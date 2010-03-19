function [OS] = change_env()
%CHANGE_ENV will handle all environmental problems
    path = pwd;
    if start_with(path, '/')
        disp 'running on LINUX system'
        addpath '~/research_linux/gnome/bioWorkspace/genomeprj/genotypeRecomb';
        addpath '~/research_linux/gnome/bioWorkspace/genomeprj/hBlockRecombo';
        addpath '~/research_linux/gnome/bioWorkspace/genomeprj/readFiles';
        addpath '~/research_linux/gnome/bioWorkspace/genomeprj/common';
        addpath '~/research_linux/gnome/bioWorkspace/genomeprj/MyDrTest';
        OS = 1;
    elseif start_with(path, 'D:')
        disp 'running on Windows, please update th path here :)';
        addpath 'D:\IUBResearch\Projects\Bioinfor\code\genotypeRecomb';
        addpath 'D:\IUBResearch\Projects\Bioinfor\code\hBlockRecombo';
        addpath 'D:\IUBResearch\Projects\Bioinfor\code\readFiles';
        addpath 'D:\IUBResearch\Projects\Bioinfor\code\common';
        addpath 'D:\IUBResearch\Projects\Bioinfor\code\MyDrTest';
        OS = 0;
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