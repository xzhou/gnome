%use mle to calculate r from genotype
function [retval] = mleR(pA, pB, n3x3)
    [m n] = size(n3x3);
    if m ~= 3 || n ~= 3
        disp 'error'
        e = MException('mleR:dim', 'in consistent dimenstion');
        throw(e);
    end
    
    pa = 1 - pA;
    pb = 1 - pB;
    
    Dmin = max(-pA*pB, -pa*pb);
    pmin = pA*pB + Dmin;
    
    Dmax = min(pA*pb, pB*pa);
    pmax = pA*pB + Dmax;
    
    %maxium likelyhood function
    %matlab support closure :-)
    function [v] = mleFunc(pAB)
        v = (2*n3x3(1,1)+n3x3(1,2)+n3x3(2,1))*log(pAB) + ...
			(2*n3x3(1,3)+n3x3(1,2)+n3x3(2,3))*log(pA-pAB) + ...
			(2*n3x3(3,1)+n3x3(2,1)+n3x3(3,2))*log(pB-pAB) + ...
			(2*n3x3(3,3)+n3x3(3,2)+n3x3(2,3))*log(1-pA-pB+pAB) + ... 
			n3x3(2,2)*log(pAB*(1-pA-pB+pAB) + (pA-pAB)*(pB-pAB));
        %NOTE: we will call fminbnb to maximize mle
        %v = -v;
    end

    [x fmin] = fminbnd(@mleFunc, Dmin, Dmax)
    
    pAB = -fmin;
    
    estD = pAB - pA*pB
    
    %D prime
    if(estD > 0)
        estDp = estD / Dmax;
    else
        estDp = estD / Dmin;
    end
    
    n = sum(sum(n3x3));
    
    corr = estD / sqrt(pA*pB*pa*pb);
    dchi = (2*n*estD^2)/(pA*pb*pB*pb);
    
    retval.r = corr;
    retval.pAB = pAB;
    retval.pA = pA;
    retval.pB = pB;
end
