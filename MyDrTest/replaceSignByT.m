function [ Rm ] = replaceSignByT( R, Threshold )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
	if nargin == 1
		Threshold = 0.1;
	end
	
	[m n] = size(R);
	
	Rm = zeros(m,n);
	
	k = 0;
	for i = 1:m
		for j = 1:n
			if abs(R(i,j)) >= Threshold
				k = k + 1;
				Rm(i,j) = R(i,j);
			else
				x = unifrnd(-1, 1);
				if x >= 0
					Rm(i,j) = abs(R(i,j));
				else
					Rm(i,j) = -abs(R(i,j));
				end
			end
		end
	end
	k
end

