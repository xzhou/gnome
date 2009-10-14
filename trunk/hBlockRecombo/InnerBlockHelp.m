classdef InnerBlockHelp
    %help functions to 
    properties
    end
    
    methods(Static)
        %% using the result to recover the sequence
        function [seq] = recoverCaseSeq(result)
            
            [nBlock repeat] = size(result);
            bestResult = cell(nBlock, 1);
            for i = 1:nBlock
                bestI = 1;
                currentBestQ = result{i,1}.finalQual;
                for j = 1:repeat
                    aResult = result{i,j};
                    if currentBestQ > aResult.finalQual
                        bestI = j;
                        currentBestQ = aResult.finalQual;
                    end
                    bestResult{i,1} = result{i,bestI};
                end
            end

            %%reconstruct
            for i = 1:nBlock
                aResult = bestResult{i,1};
                aSeq = blockReconstruct(aResult.refBlock, aResult.finalFreq);   
                if i == 1
                    seq = aSeq;
                else
                    seq = [seq aSeq];
                end
            end
        end
        function [signRate] = calcSignRate(caseR, refR)
            [mc nc] = size(caseR);
            [mr nr] = size(refR);
            if mc ~= mr || nc ~= nr
                e = MException('InnerBlockHelp:calcSingRate', 'inconsistent');
                throw(e);
            end
            
            caseSign = sign(caseR);
            refSign = sign(refR);
            
            mix = caseSign.*refSign;
            mix(logical(eye(mc))) = -2;
            n = sum(sum(mix>=0));
            signRate = n*1.0/mc/(mc-1);
        end
    end
end

