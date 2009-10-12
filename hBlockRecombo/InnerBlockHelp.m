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
    end
end

