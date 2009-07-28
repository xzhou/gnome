function [ E ] = evaluatePop( learnedSample, targetRS, targetF, K)
%evaluatePop calculate the evaluation value of a haplotype
% @learnedSample: the haplotype of the a learned Sample
% @targetR: the target R square values
% @targetF: the target single allele frequencies
% @K: is a constant of rs diff and the single allele freq diff

	if nargin == 3
		K = 0.1;	%the default value of K is 0.1

	learnedSample = majorize(learnedSample);
	
	[r, rs] = calculateR(learnedSample);
	F = calculateSinlgeAlleleFreq(learnedSample);
	
	rsDiff = (rs - targetRS).*(rs - targetRS);
	
	rsError = nansum(rsDiff(:))/2;
	
	FDiff = (F - targetF).*(F - targetF);
	
	FError = nansum(FDiff(:))/2;
	
	E = (1-K)*rsError + K*FError;
		
end

