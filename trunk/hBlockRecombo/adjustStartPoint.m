function newCurrentSeq = adjustStartPoint(currentSeq, refSeq, blocks)
    nS = length(refSeq(:,1));
    Len = length(refSeq(1,:));
    newCurrentSeq = zeros(nS, Len);
    [nBlock tmp] = size(blocks);
    currentBlockFreqInfo = cell(nBlock, 1);
    refBlockFreqInfo = cell(nBlock, 1);
    parfor i = 1:nBlock
        refBlockFreqInfo{i,1} = getBlockFreq(refSeq, blocks(i,:));
    end
    parfor i = 1:nBlock
        currentBlockFreqInfo{i,1} = getBlockFreq(currentSeq, blocks(i,:));
    end
    currentMatchedRef = blockCheck(currentBlockFreqInfo, refBlockFreqInfo, blocks);
    %currentMatchedCase = blockCheck(currentBlockFreqInfo,
    %caseBlockFreqInfo, blocks);
    
    for i = 1:nBlock
        currentMatchedRef{i, 1} = currentMatchedRef{i,1}(1:end-2, :);
        %for each
        for j = 1 : nS
            tempRefSeq4 = refSeq(:, blocks(i,1):blocks(i,2));
            k = 1; 
            found = 0;
            while ((k<=length(currentMatchedRef{i,1}(:,1))&&(found == 0)))
                found = isequal(tempRefSeq4(j,:), currentMatchedRef{i,1}(k, 1:end-2));
                k = k+1;
            end
            if ((found == 1)&&(currentMatchedRef{i,1}(k-1, end-1)~=0))
                newCurrentSeq(j, blocks(i,1):blocks(i,2)) = currentMatchedRef{i,1}(k-1, 1:end-2);
                currentMatchedRef{i,1}(k-1, end-1) = currentMatchedRef{i,1}(k-1, end-1)-1;
                currentMatchedRef{i,1}(k-1, end) = currentMatchedRef{i,1}(k-1, end)-1;
            else
                disBetweenSeq = zeros(length(currentMatchedRef{i,1}), 1);
                for l = 1 : length(currentMatchedRef{i,1}(:,1))
                    if ((currentMatchedRef{i,1}(l,end-1))>(currentMatchedRef{i,1}(l, end)))
                        disBetweenSeq(l,1) = sum(tempRefSeq4(j,:) == currentMatchedRef{i,1}(l, 1:end-2));
                    end
                end
                [maxDisSum,maxDisIdx] = max(disBetweenSeq(:,1));
                newCurrentSeq(j, blocks(i,1):blocks(i,2)) = currentMatchedRef{i,1}(maxDisIdx, 1:end-2);
                currentMatchedRef{i,1}(maxDisIdx, end-1) = currentMatchedRef{i,1}(maxDisIdx, end-1)-1;
            end
        end
    end
    
end