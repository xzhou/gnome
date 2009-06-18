# TODO: Add comment
# 
# Author: xzhou
###############################################################################

#function calcAllR will calculate the r value for all pair of genotypes 
calcAllR <- function(g1,...)
{
	gvars <- sapply( g1, function(x) (is.genotype(x) && nallele(x)==2) )
	if(any(gvars==FALSE))
	{
		warning("Non-genotype variables or genotype variables ",
				"with more or less than two alleles detected. ",
				"These variables will be omitted: ",                
				paste( colnames(g1)[!gvars] , collapse=", " )
		)
		g1 <- g1[,gvars]
	}
	
	
	
	#do some thing here, I must try to recover those error
	#those column is invalid
	fnames = which(!gvars)
	
	P <- matrix(nrow=ncol(g1),ncol=ncol(g1))
	rownames(P) <- colnames(g1)
	colnames(P) <- colnames(g1)
	
	
	#add by xzhou@indiana, dump pA, pB, pAB
	P <- D <- Dprime <- nobs <- chisq <- p.value <- corr <- R.2 <- pA <- nameA <- namea <- pB <- nameB <- nameb <- pAB <- cAMax <-cAMin <- cANA <-cBMax <-cBMin <- cBNA <- P
	P <- n11 <- n12 <- n13 <- n21 <- n22 <- n23 <- n31 <- n32 <- n33 <- P
	P <- R.1 <- P
	
	cat("@LD(g1...) ncol=", ncol(g1))
	cat("@LD(pA...) ncol=", ncol(pA))
	
	for(i in 1:(ncol(g1)-1) )
		for(j in (i+1):ncol(g1) )
		{
			ld <- LD( g1[,i], g1[,j] )
			
			D      [i,j] <- ld$D
			Dprime [i,j] <- ld$"D'"
			
			
			if(FALSE && abs(ld$"r") <= 0.1)
			{
				corr   [i,j] <- 1.0
			}
			else
			{
				corr   [i,j] <- ld$"r"
			}
			
			#cut lower triangle
			if(FALSE)
			{
				corr   [i,j] <- ld$"r"
				if((j-i)>=40)
				{
					corr   [i,j] = 1.0
				}
			}
			
			R.2    [i,j] <- ld$"R^2"
			nobs   [i,j] <- ld$"n"
			chisq  [i,j] <- ld$"X^2"
			p.value[i,j] <- ld$"P-value"
			pA[i,j] <- ld$"pA"
			cAMax[i, j] <- ld$"cAMax"
			cAMin[i, j] <- ld$"cAMin"
			cANA[i, j] <- ld$"cANA"
			
			pB[i,j] <- ld$"pB"
			cBMax[i,j] <- ld$"cBMax"
			cBMin[i,j] <- ld$"cBMin"
			cBNA[i,j] <- ld$"cBNA"
			
			pAB[i,j] <- ld$"pAB"
			
			nameA[i,j] <- ld$"mA"
			nameB[i,j] <- ld$"mB"
			namea[i,j] <- ld$"ma"
			nameb[i,j] <- ld$"mb"
			
			n11[i,j] <- ld$"n11"
			n12[i,j] <- ld$"n12"
			n13[i,j] <- ld$"n13"
			n21[i,j] <- ld$"n21"
			n22[i,j] <- ld$"n22"
			n23[i,j] <- ld$"n23"
			n31[i,j] <- ld$"n31"
			n32[i,j] <- ld$"n32"
			n33[i,j] <- ld$"n33"
		}
	
	#out put the value
	
	
	retval <- list(
			call=match.call(),
			"D"=D,
			"D'"=Dprime,
			"r" = corr,
			"R^2" = R.2,
			"n"=nobs,
			"X^2"=chisq,
			"P-value"=p.value,
			"pA"=pA,
			"pB"=pB,
			"pAB"=pAB,
			"cAMax"=cAMax,
			"cAMin"=cAMin,
			"cANA"=cANA,
			"cBMax"=cBMax,
			"cBMin"=cBMin,
			"cBNA"=cBNA,
			"nameA"=nameA,
			"nameB"=nameB,
			"namea"=namea,
			"nameb"=nameb,
			
			"n11" = n11,
			"n12" = n12,
			"n13" = n13,
			"n21" = n21,
			"n22" = n22,
			"n23" = n23,
			"n31" = n31,
			"n32" = n32,
			"n33" = n33,
			"fnames" = fnames
	)
	
	class(retval) <- "LD.data.frame"
	
	retval
}
	
	
rcalc <- function(g1, g2, ...)
{
	if(nallele(g1) != 2 || nallele(g2) != 2)
	{
		#warning("Program only support 2-allele genotype")
		#g1
		#g2
		retval <- list("r" = 0.0, "pAB" = 0.0, "pA" = 0.0, "pB" = 0.0)
		return(retval)
	}
	
	prop.A <- summary(g1)$allele.freq[,2]
	prop.B <- summary(g2)$allele.freq[,2]
	
	count.A <- summary(g1)$allele.freq[,1]
	count.B <- summary(g2)$allele.freq[,1]
	
	
	major.A <- names(prop.A)[which.max(prop.A)]
	minor.A <- names(prop.A)[which.min(prop.A)]
	major.B <- names(prop.B)[which.max(prop.B)]
	minor.B <- names(prop.B)[which.min(prop.B)]
	
	pA <- max(prop.A, na.rm=TRUE)
	cAMax <- count.A[major.A]
	cAMin <- count.A[minor.A]
	cANA <- count.A["NA"]
	
	pB <- max(prop.B, na.rm=TRUE)
	cBMax <- count.B[major.B]
	cBMin <- count.B[minor.B]
	cBNA <- count.B["NA"]
	
	pa <- 1-pA
	pb <- 1-pB
	
	Dmin <- max(-pA*pB, -pa*pb)
	pmin <- pA*pB + Dmin;
	
	Dmax <- min(pA*pb, pB*pa);
	pmax <- pA*pB + Dmax;
	
	counts <- table(
			allele.count(g1, major.A),
			allele.count(g2, major.B)
	)
	
	#print(counts)
	
	n3x3 <- matrix(0, nrow=3, ncol=3)
	colnames(n3x3) <- rownames(n3x3) <- 0:2
	
	# ensure the matrix is 3x3, with highest frequency values in upper left
	for(i in rownames(counts))
		for(j in colnames(counts))
			n3x3[3-as.numeric(i),3-as.numeric(j)] <- counts[i,j]
	
	
	loglik <- function(pAB,...)
	{
		(2*n3x3[1,1]+n3x3[1,2]+n3x3[2,1])*log(pAB) +
				(2*n3x3[1,3]+n3x3[1,2]+n3x3[2,3])*log(pA-pAB) +
				(2*n3x3[3,1]+n3x3[2,1]+n3x3[3,2])*log(pB-pAB) +
				(2*n3x3[3,3]+n3x3[3,2]+n3x3[2,3])*log(1-pA-pB+pAB) + 
				n3x3[2,2]*log(pAB*(1-pA-pB+pAB) + (pA-pAB)*(pB-pAB))
	}
	
	# SAS code uses:
	#
	#s <- seq(pmin+0.0001,pmax-0.0001,by=0.0001)
	#lldmx <- loglik(s)
	#maxi <- which.max(lldmx)
	#pAB <- s[maxi]
	
	# but this should be faster:
	solution <- optimize(
			loglik,
			lower=pmin+.Machine$double.eps,
			upper=pmax-.Machine$double.eps,
			maximum=TRUE
	)
	pAB <- solution$maximum
	
	estD <- pAB - pA*pB
	if (estD>0)  
		estDp <- estD / Dmax
	else
		estDp <- estD / Dmin
	
	n <-  sum(n3x3)
	
	corr <- estD / sqrt( pA * pB * pa * pb )
	
	dchi <- (2*n*estD^2)/(pA * pa * pB* pb)
	dpval <- 1 - pchisq(dchi,1)
	
	#return r
	
	retval <- list("r" = corr, "pAB" = pAB, "pA" = pA, "pB" = pB)
	
	retval
}
