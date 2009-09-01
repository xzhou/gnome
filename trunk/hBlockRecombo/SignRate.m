function [rate totalSign correctSign] = SignRate(targetR, sampleR, blockLens)
%calculate the sign recover rate for the blocks
% @targetR the target R value
% @sampleR the sample R value
%
    if nargin == 2
        targetSign = sign(targetR);
        sampleSign = sign(sampleR);

        signDiff = targetSign.*sampleSign;
        signDiff(logical(eye(size(signDiff)))) = 0;

        correctSign = sum(sum(double(signDiff == 1)))/2;

        [m n] = size(targetR);

        totalSign = m*(m-1)/2;

        rate = correctSign/totalSign;
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