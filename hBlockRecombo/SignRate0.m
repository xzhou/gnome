function [rate totalSign correctSign] = SignRate0(targetR, sampleR, blockLens)
%calculate the sign recover rate for the blocks
% @targetR the target R value
% @sampleR the sample R value
%
    if nargin == 2
        targetR(logical(targetR==0))=NaN;
        sampleR(logical(sampleR==0))=NaN;
        targetSign = sign(targetR);
        sampleSign = sign(sampleR);
        signDiff = targetSign.*sampleSign;
        signDiff(logical(eye(size(signDiff)))) = 0;
        correctSign = nansum(nansum(double(signDiff == 1)))/2;
        [m ~] = size(targetR);
        testR=targetR;
        testR(logical(triu(ones(m))))=1;
        adjust = sum(sum(testR==0));
        totalSign = m*(m-1)/2;
        %remove Nan value
        totalSign = totalSign - adjust;
        rate = (correctSign-adjust)/totalSign;
    else
        nBlock = length(blockLens);
        maskMatrix = zeros(size(targetR));
        %remove the sign inside block
        for i = 1:nBlock
            if i == 1
                lim1 = 1;
                lim2 = blockLens(i);
            else
                lim1 = sum(blockLens(1:i-1));
                lim2 = sum(blockLens(1:i));
            end
            maskMatrix(lim1:lim2, lim1:lim2) = 1;
        end
        maskMatrix = logical(maskMatrix);
        reverseMask = ~maskMatrix;
        targetSign = sign(targetR);
        sampleSign = sign(sampleR);
        
        signDiff = targetSign.*sampleSign;
        signDiff(maskMatrix) = 0;
        
        correctSign = sum(sum(double(signDiff == 1)))/2;
        totalSign = sum(sum(double(reverseMask)))/2;
        rate = correctSign/totalSign;
    end
end