function [compareResult] = commonGenotype(commonStart, commonEnd)
    [m n] = size(commonStart);
    compareResult = cell(3,1);
    for i = 1:m
        aResult = [];
        sBlock = commonStart{i,1};
        eBlock = commonEnd{i,1};
        [a b] = size(sBlock);
        for j = 1:a
            k = findBlock(sBlock(j,2), eBlock);
            if k > 0
                aResult = [aResult; [sBlock(j,:), eBlock(k,:)]];
            end
        end
        compareResult{i, 1} = aResult;
    end
end

function [index] = findBlock(x, aList)
    index = 0;
    [m n] = size(aList);
    for i = 1:m
        if aList(i,2) == x
            index = i;
            return;
        end
    end
    return;
end
