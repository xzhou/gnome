function indxtomiss2 = indxshift(indxtomiss2, indxtomiss1)
	% in the original array, both indxtomiss1 and indxtomiss2 are to be removed
	% if we first remove indxtomiss1, then the indxtomiss2 should be modified
	% in order to operate on the new array
	
	[x1,y1] = size(indxtomiss2);
	[x2,y2] = size(indxtomiss1);
	
	if x1 == 0 && x2 == 0 && y1 == 0 && y2 == 0
		indxtomiss2 = [];
		return
	end
	
	indxtomiss2 = setdiff(indxtomiss2, indxtomiss1); %sorting and remove redundency
	[m,n]=size(indxtomiss1);
	if n == 1
		indxtomiss1 = indxtomiss1';
	end
	[m,n]=size(indxtomiss2);
	if n == 1
		indxtomiss2 = indxtomiss2';
	end
	Len2 = length(indxtomiss2);
	[tmp indx] = sort([indxtomiss2 indxtomiss1]);
	indxtomiss2 = indxtomiss2 - (find(indx<=Len2) - [1:Len2]);