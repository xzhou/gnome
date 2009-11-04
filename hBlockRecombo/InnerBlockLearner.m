classdef InnerBlockLearner
    %class InnerBlockLearningResult record the learning result
    properties
        
    end
    
    methods
        function [seq] = innerBlockDriver(caseSeq4, refSeq4, blocks)    
        %cd 'D:\IUBResearch\Projects\Bioinfor\data\88_77_CEU_YRI_DATA';
            [nBlock tmp] = size(blocks);

            nS = length(caseSeq4);
            Len = length(caseSeq4(1,:));

            caseBlockFreqInfo = cell(nBlock, 1);

            parfor i = 1:nBlock
                caseBlockFreqInfo{i,1} = getBlockFreq(caseSeq4, blocks(i,:));
            end

            refBlockFreqInfo = cell(nBlock, 1);
            parfor i = 1:nBlock
                refBlockFreqInfo{i,1} = getBlockFreq(refSeq4, blocks(i,:));
            end

            matchedCase = blockCheck(caseBlockFreqInfo, refBlockFreqInfo, blocks);
            refMatchedCase = blockCheck(refBlockFreqInfo, caseBlockFreqInfo, blocks);

            startParallel(2);

            nRepeat = 10;
            result = cell(nBlock,nRepeat);

            for i=1:nBlock
                block = blocks(i,:);
                a = block(1, 1);
                b = block(1, 2);
                refBlock = refMatchedCase{i,1}(1:end-2,1:end-2);
                refFreq = refMatchedCase{i,1}(1:end-2, end-1);
                refCaseFreq = refMatchedCase{i,1}(1:end-2, end);
                caseBlock = matchedCase{i,1}(1:end-2, 1:end-2);
                caseFreq = matchedCase{i,1}(1:end-2, end-1);
                caseRefFreq = matchedCase{i,1}(1:end-1, end);
                refSeq = refSeq4(:,a:b);
                caseSeq = caseSeq4(:,a:b);
                alleleMapping = getMajorAllele(refSeq);
                fprintf(1, '\n************ learning block %d ***********\n', i);
                parfor k = 1:nRepeat
                    [aResult] = innerBlockLearning(caseBlock, caseFreq, refBlock, refFreq, refCaseFreq, alleleMapping);
                    [bResult] = innerBlockLearning(caseBlock, caseFreq, refBlock, refFreq, refCaseFreq, alleleMapping,1);
                    try
                        result{i, k} = aResult;
                    catch exception
                        fprintf(1, 'error\n');
                        rethrow(exception);
                    end
                    fprintf(1, 'block = %d repeat = %d \t a = %f,\t b = %f\n',i, k, aResult.fDistance, bResult.fDistance);
                end
            end
            %% reconstruct the learned block we use the
            seq = InnerBlockHelp.recoverCaseSeq(result);
            
            %save everything for post analysis
            save('innterseq.mat', 'seq');
            save('innerBlockLearner.mat');
        end
        function [result] = innerBlockLearning(caseBlock, caseFreq, refBlock, refFreq, refcaseFreq, alleleMapping, alpha)
            %alpha is the weight of r square in the learing process
            disp 'hello';
            if nargin == 5
                alleleMapping = getMajorAllele(refSeq4);
                alpha = 0.6;
            elseif nargin == 6
                alpha = 0.6;
            end

            %defines the randomness of searching strategy
            expT = 1.0e-5;

            %% learning algorithm
            maxIt = 1e0;
            itr = 0;

            targetSeq = blockReconstruct(caseBlock, caseFreq);
            caseR = calcR(targetSeq, alleleMapping);
            caseRs = caseR.*caseR;
            caseAlleleFreq = GnomeCalculator.getSingleAlleleFreq(targetSeq, alleleMapping);

            currentSeq = blockReconstruct(refBlock, refFreq);
            refR = calcR(currentSeq, alleleMapping);
            currentR = refR;
            currentRs = currentR.*currentR;
            currentAlleleFreq = GnomeCalculator.getSingleAlleleFreq(currentSeq, alleleMapping);

            %initialize
            currentFreq = refFreq;
            currentQuality = eval(caseRs, currentR.*currentR, caseAlleleFreq, currentAlleleFreq, alpha);

            %fprintf(1, '%f\n', currentQuality);

            historySize = 1000;
            previousKQuality = zeros(historySize, 1);

            while itr < maxIt
                itr = itr + 1;
                [newSeq4 newFreq] = getNextSeq(refBlock, currentFreq);
                newR = calcR(newSeq4, alleleMapping);
                newRs = newR.*newR;
                %newP is the single allele frequency of the new allele frequency
                newP = GnomeCalculator.getSingleAlleleFreq(newSeq4, alleleMapping);
                newQuality = eval(caseRs, newRs, caseAlleleFreq, newP, alpha);
                Qdiff = newQuality - currentQuality;
                %using statistic hill climbing algorithm to change
                if Qdiff < 0
                    p = 1/(1+exp(Qdiff/expT));
                else
                    p = 1/(1+exp(Qdiff/expT));
                end
                x = rand(1);
                if x < p
                    currentFreq = newFreq;
                    currentQuality = newQuality;
                    currentRs = newRs;
                    currentR = newR;
                    currentP = newP;
                    %fprintf(1, 'itr = %d\t newQuality = %.15f\n', itr, newQuality);
                else
                    %do nothing
                end
            end

            %end of learning
            result.finalFreq = currentFreq;
            result.finalQual = currentQuality;
            result.finalR = currentR;
            result.refFreq = refFreq;
            result.refcaseFreq = refcaseFreq;
            result.refBlock = refBlock;

            %calculate the frequency distance
            a = sum(abs(result.refcaseFreq - result.finalFreq));

            b = sum(abs(result.refcaseFreq - result.refFreq));

            %if fDistance < 1, the learning process is good
            result.fDistance = a*1.0/b;
        end
        function [newSeq, newBlockFreq] = getNextSeq(hyplotypes, currentBlockFreq)
            newBlockFreq = blockNaiveMutate(currentBlockFreq);
            newSeq = blockReconstruct(hyplotypes, newBlockFreq);
        end
        %evaluate normalized r square difference and p difference 
        function [newQuality] = eval(targetRs, currentRs, targetFreq, currentFreq, alpha)
            %we think R square has more weight than single allele frequence
            if nargin == 4
                alpha = 0.6;
            end

            pDiff = abs(targetFreq - currentFreq);
            normalDiff = sum(pDiff)/length(targetFreq);   %normalize

            rDiff = calcRDiff(targetRs, currentRs);
            [m, n] = size(targetRs);
            if m ~= n
                e = MException('blockReconstruct:x', 'in consistent dimenstion');
                throw(e);
            end
            nElements = m*(m-1)/2;
            normalRDiff = rDiff/nElements;

            newQuality = alpha*normalRDiff + (1-alpha)*normalDiff;
            newQuality = 100*newQuality;
        end
        function [newQuality] = calcRDiff(targetRs, newRs, smallRFilter)

            if nargin == 2
                smallRFilter = 0.0;
            end
            %make sure remove the diagnal elements
            targetRs(logical(eye(size(targetRs)))) = 0;
            newRs(logical(eye(size(newRs)))) = 0;

            %filter the small r values
            targetRs(targetRs < smallRFilter) = 0;
            newRs(targetRs < smallRFilter) = 0;

            newQuality = sum(sum(abs(targetRs - newRs)))/2;
        end
    end
end

