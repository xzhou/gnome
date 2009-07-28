function [ c0x ] = calculateSinlgeAlleleFreq( haplotype )
%CALCULATESINGLEALLELEFREQ calculate the single allele frequency
% haplotype is the phased data
	haplotype = majorize(haplotype);
	
	c0x = sum(haplotype == 0);
	
end

