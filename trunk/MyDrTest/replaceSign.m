function [ rValue ] = replaceSign( R, recoveredSigns)
%REPLACESIGN replace the single r and double r
%   Detailed explanation goes here
%@Param F is the single allele frequency
%@Param R is the pairwise frequency
%@Recovered Signs is a 0 1 -1 matrix indicating the sign, 0 as 
% unknown signs, 1 as "+" and -1 as "-"
	
	%rValue = abs(R).*recoveredSigns
	
	[m1, n1] = size(R);
	[m2, n2] = size(recoveredSigns);
	
	%terminate the current program
	if m1 ~= m2 || n1 ~= n2
		exit;
	end
	
	for i = 1:m1
		for j = 1:n1
			%if the sign is unknown, use random sign
			if recoveredSign(i,j) == 0
				x = unifrnd(-1, 1)
				if x >= 0
					rValue(i,j) = R(i,j);
				else
					rValue(i,j) = R(i,j) * -1;
				end
			elseif recoveredSign(i,j) == 1
				rValue(i,j) = abs(R(i,j));
			elseif recoveredSign(i,j) == -1
				rValue(i,j) = abs(R(i,j))*-1
			end
		end
	end
	fValue = F;
end

