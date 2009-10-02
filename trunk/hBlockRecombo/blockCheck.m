function [ caseSeq4WithFreq ] = blockCheck( caseSeq4WithFreq, refSeq4WithFreq, block)
%function blockCheck will check what hyplotype block are in reference and
%return this information to the last column
	[mBlock, tmp] = size(block);
	
	for i = 1:mBlock
		caseBlock = caseSeq4WithFreq{i,1};
		refBlock = refSeq4WithFreq{i,1};
		
		[nHyplotype, len] = size(caseBlock);
		
		blockMatches = zeros(nHyplotype, 1);
		
		%m is the number of different case hyplotype
		
		for k = 1:nHyplotype
			aHyplotype = caseBlock(k, :);
			refFreq = findHyplotype(refBlock, aHyplotype);
			blockMatches(k,1) = refFreq;
		end
		
		caseBlock = [caseBlock blockMatches];
		caseSeq4WithFreq{i,1} = caseBlock;
	end
	caseSeq4WithFreq
end

function [refFreq] = findHyplotype(refSeq4WithFreq, patternHyplotypeWithFreq)
	[m n] = size(refSeq4WithFreq);
	refFreq = 0;
	for i = 1:m
		aHyplotype = refSeq4WithFreq(i,:);
		if sum(patternHyplotypeWithFreq(:,1:end-1) ~= aHyplotype(:,1:end-1)) == 0
			patternHyplotypeWithFreq
			aHyplotype
			refFreq = aHyplotype(1, end);
			break;
		end
	end
end

