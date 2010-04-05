function [above95 above99] = getRealPersonNo(nS, StatS, SortStatS)

	Matrix95 = StatS>SortStatS(int16(nS*0.95));
	Matrix99 = StatS>SortStatS(int16(nS*0.99));

	above95 = sum(StatS>SortStatS(int16(nS*0.95)));
	above99 = sum(StatS>SortStatS(int16(nS*0.99)));
	
	for i=1:2:nS
		if((Matrix95(i,1).*Matrix95(i+1,1))==1)
			above95=above95-1;
		end
		if((Matrix99(i,1).*Matrix99(i+1,1))==1)
			above99=above99-1;
		end
		
	end
end
	