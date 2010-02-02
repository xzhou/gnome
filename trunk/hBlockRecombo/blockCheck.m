function [ caseSeq4WithFreq ] = blockCheck(caseSeq4WithFreq, refSeq4WithFreq, block)
%function blockCheck will check what hyplotype block are in reference and
%return this information to the last column

    if nargin == 2
        [a b] = size(caseSeq4WithFreq)
        block = [1,b];
    end
    
	[mBlock, tmp] = size(block); 
    
	for i = 1:mBlock
		caseBlock = caseSeq4WithFreq{i,1};
		refBlock = refSeq4WithFreq{i,1};
        
        %get block frequency
		alleleMapping = getMajorAllele(refBlock(:,1:end-1));
        caseFreq = getSingleAlleleFreq(caseBlock(:,1:end-1), alleleMapping);
        caseFreq = [caseFreq, 0, 0];

		[nHyplotype, len] = size(caseBlock);
		
		blockMatches = zeros(nHyplotype, 1);
		
		%m is the number of different case hyplotype
		
		for k = 1:nHyplotype
			aHyplotype = caseBlock(k, :);
			refFreq = findHyplotype(refBlock, aHyplotype);
			blockMatches(k,1) = refFreq;
		end
		
		caseBlock = [caseBlock blockMatches];
        caseBlock = [caseBlock;caseFreq];
        caseBlock = [caseBlock; [alleleMapping, 0, 0]];
        
		caseSeq4WithFreq{i,1} = caseBlock;
	end
	caseSeq4WithFreq;

end

function [singleAlleleFreq] = getSingleAlleleFreq(seq4, majorAlleleMapping)
    if nargin == 1
        majorAlleleMapping = getMajorAllele(seq4);
    end
    
    [nS, tmp] = size(seq4);
    
    int2S = (seq4 == repmat(majorAlleleMapping,nS,1)) + 0;
    singleAlleleFreq = sum(int2S);
end

function [refFreq] = findHyplotype(refSeq4WithFreq, patternHyplotypeWithFreq)
	[m n] = size(refSeq4WithFreq);
	refFreq = 0;
	for i = 1:m
		aHyplotype = refSeq4WithFreq(i,:);
		if sum(patternHyplotypeWithFreq(:,1:end-1) ~= aHyplotype(:,1:end-1)) == 0
			patternHyplotypeWithFreq;
			aHyplotype;
			refFreq = aHyplotype(1, end);
			break;
		end
	end
end

