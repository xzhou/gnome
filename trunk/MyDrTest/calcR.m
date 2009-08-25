function [r p] = calcR(pedFileName)
    cd '/home/xzhou/research_linux/gnome/workspace/data/realdata'
    
    caseData = importdata('caseonly.ped');
    controlData = importdata('controlonly.ped');
    
    seqS = caseData(6:end,:);
    seqR = controlData(6:end,:);
    
    
    
end