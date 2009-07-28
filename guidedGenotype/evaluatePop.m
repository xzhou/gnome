function [ E ] = evaluatePop( learnedSample, targetRS, targetF )
%evaluatePop calculate the evaluation value of a haplotype
% @learnedSample: the haplotype of the a learned Sample
% @targetR: the target R square values
% @targetF: the target single allele frequencies

	[r, rs] = calculateR(learnedSample);
	
	

end

