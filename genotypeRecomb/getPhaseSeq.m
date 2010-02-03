function [phasedSeq] = getPhaseSeq(id, fastaSeq)
    [fm fn] = size(fastaSeq);
    [im in] = size(id);
    
    phasedSeq = [];
    for i = 1:im
        idName = id{i};
        flag = 0;
        for j = 1:fm
            fastaIdName = fastaSeq(j,1).Header;
            k = findstr(idName, fastaIdName);
            if ~isempty(k)
                phasedSeq = [phasedSeq; fastaSeq(j,1); fastaSeq(j+1,1)];
                flag = 1;
                break;
            end
        end
        if flag == 0
            fprintf(1, '%s can not find\n', idName);
        end
    end
end