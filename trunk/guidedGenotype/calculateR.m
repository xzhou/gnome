function [ rValues, pA ] = calculateR( haplotype, func)

	if nargin == 1, fun=@rsquare; end  % fun in default is @rsquare

	[row,col]=size(A);

	for m=1:(col-1)
		for n=(m+1):col
			x11=0;
			x12=0;
			x21=0;
			x22=0;
			for i=1:row
				if haplotype(i,m)=='0' && haplotype(i,n)=='0'
					x11=x11+1;
				elseif haplotype(i,m)=='0' && haplotype(i,n)=='1'
					x12=x12+1;
				elseif haplotype(i,m)=='1' && haplotype(i,n)=='0'
					x21=x21+1;
				elseif haplotype(i,m)=='1' && haplotype(i,n)=='1'
					x22=x22+1;
				end
			end
			[rValue, rsValue] = feval(fun,x11,x12,x21,x22);
			rs(m,n) = rsValue
			rs(n,m)=rs(m,n);

			r(m,n) = rValue;
			r(n,m) = rValue;
		end
	end
	[row,col]=size(rs);

	for i=1:row
		rs(i,i)=sum(rs(:,i))/(row-1); % rs(i,i) has the averaged LD for ith SNP
		r(i,i) = NaN;
	end
	
end

