function [ haplotype ] = majorize( haplotype )
%MAJORIZE change the encoding so that 0 is the major and 1 is the minor
%   the unmajorized 01 encoded phased data
	[nrow, ncol] = size(haplotype)
	
	
	for j = 1:ncol
		aSnp = haplotype(:,j);
		count0 = sum(aSnp == 0);
		count1 = sum(aSnp == 1);
		
		if count0 < count1
			%change the major and minor
			aSnp(aSnp == 0) = -1;	%save the 0s
			aSnp(aSnp == 1) = 0;	%flip
			aSnp(aSnp == -1) = 1;	%flip
			haplotype(:,j) = aSnp;
		end
	end
end

