function [common] = compareCaseRef(caseInfo, refInfo, blocks)
    [m n] = size(blocks);
    common = cell(3,1);
    for i = 1:m
        aCaseInfo = caseInfo{i, 1};
        aRefInfo = refInfo{i, 1};
        caseBlockType = aCaseInfo.uniqueBlock;
        refBlockType = aRefInfo.uniqueBlock;
        a = [];
        [nType, tmp] = size(caseBlockType);
        for j = 1:nType
            aTypeSeq = caseBlockType(j,:);
            x = containsSeq(aTypeSeq, refBlockType);
            if x > 0
                a = [a; j, x, caseInfo{i,1}.freq(j), refInfo{i,1}.freq(x)];
            end
        end
        common{i,1} = a;
    end
end

function [t] = containsSeq(aSeq, seqs)
    [m n] = size(seqs);
    for i = 1:m
        if sum(aSeq ~= seqs(i,:)) == 0
            t = i;
            return;
        end
    end
    t = 0;
end