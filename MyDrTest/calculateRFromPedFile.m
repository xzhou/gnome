function [r p] = calculateRFromPedFile(pedFileName, ncol)
    
    if nargin == 0
        cd '/home/xzhou/research_linux/gnome/workspace/data/realdata'
        pedFileName = '80SNPHapmap.ped'
        ncol = 166;
    end

    %try
        fid = fopen(pedFileName);
        allInfo = textscan(fid, ['%s' repmat(' %s', 1, ncol-1)]);
        genodata = allInfo(:,7:end);
        
        nrow = length(genodata{:,1});
        ncol = length(genodata);
        
        %convert to sequence
        allSeq = {};
        for i = 1:nrow
            aSeq = [];
            for j = 1:ncol
                aSeq = strcat(aSeq, genodata{1,j}(i,1));
            end
            aSeq = aSeq{1};
            allSeq{i,1} = aSeq;
        end
        
        nrow = length(allSeq);
       
        int4R = zeros(nrow, ncol);
        for i = 1:nrow
            int4R(i,:) = nt2int(allSeq{i,1}) - 1;
        end
        

end