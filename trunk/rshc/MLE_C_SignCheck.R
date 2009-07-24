# Check if the sign of MLE recovered C has a high sign recover rate
# 
# Author: xzhou
###############################################################################

source("SHC.R")
library(mygenetics)
# read the dump R value
readDumpRValue <- function()
{
	x <- read.table(file = "../data/sim_4000seq/80SNP_CEU_sim_4000seq.freq.txt", 
			header = TRUE)
	
	c = x[, c(1, 2, 18)]
	
	nSnps = max(c[,2])
	
	
	m = ncol(c)
	
	dumpR = matrix(0.0, nSnps, nSnps)
	for(k in 1:m)
	{
		i = c[k, 1]
		j = c[k, 2]
		r = c[k, 3]
		dumpR[i,j] <- dumpR[j,i] <- r
	}
	diag(dumpR) <- NA
	dumpR
}

readMLECounts <- function(verbose = F)
{
	x = read.table("../convexAnalysis/c11m.out")
	m = nrow(x)
	n = ncol(x)
	c11m = matrix(as.numeric(as.matrix(x)), m, n)
	if(verbose)
		print(c11m)
	
	x = read.table("../convexAnalysis/c12m.out")
	m = ncol(x)
	n = ncol(x)
	c12m = matrix(as.numeric(as.matrix(x)), m, n)
	
	x = read.table("../convexAnalysis/c13m.out")
	m = ncol(x)
	n = ncol(x)
	c13m = matrix(as.numeric(as.matrix(x)), m, n)
	
	x = read.table("../convexAnalysis/c21m.out")
	c21m = matrix(as.numeric(as.matrix(x)), m, n)
	
	x = read.table("../convexAnalysis/c22m.out")
	m = ncol(x)
	n = ncol(x)
	c22m = matrix(as.numeric(as.matrix(x)), m, n)
	
	x = read.table("../convexAnalysis/c23m.out")
	m = ncol(x)
	n = ncol(x)
	c23m = matrix(as.numeric(as.matrix(x)), m, n)
	
	x = read.table("../convexAnalysis/c31m.out")
	m = ncol(x)
	n = ncol(x)
	c31m = matrix(as.numeric(as.matrix(x)), m, n)
	
	x = read.table("../convexAnalysis/c32m.out")
	m = ncol(x)
	n = ncol(x)
	c32m = matrix(as.numeric(as.matrix(x)), m, n)
	
	x = read.table("../convexAnalysis/c33m.out")
	m = ncol(x)
	n = ncol(x)
	c33m = matrix(as.numeric(as.matrix(x)), m, n)
	
	pA = read.table("../convexAnalysis/pA")
	
	ret = list(
			"c11m" = c11m, 
			"c12m" = c12m,
			"c13m" = c13m,
			"c21m" = c21m,
			"c22m" = c22m,
			"c23m" = c23m,
			"c31m" = c31m,
			"c32m" = c32m,
			"c33m" = c33m, 
			"pA" = pA
					)
	ret
}


#' using counts to calculates the R value
#' tested correct
calculateRValue <- function(pA, pB, c11m, c12m, c13m, c21m, c22m, c23m, c31m, c32m, c33m)
{
	pa <- 1 - pA
	pb <- 1 - pB
	
	Dmin <- max(-pA*pB, -pa*pb)
	pmin <- pA*pB + Dmin;
	
	Dmax <- min(pA*pb, pB*pa);
	pmax <- pA*pB + Dmax;
	
	#genotype counts
	n3x3 <- matrix(0, 3, 3)
	
	n3x3[1, 1] = c11m
	n3x3[1, 2] = c12m
	n3x3[1, 3] = c13m
	n3x3[2, 1] = c21m
	n3x3[2, 2] = c22m
	n3x3[2, 3] = c23m
	n3x3[3, 1] = c31m
	n3x3[3, 2] = c32m
	n3x3[3, 3] = c33m
	
	loglik <- function(pAB, ...)
	{
		(2*n3x3[1,1]+n3x3[1,2]+n3x3[2,1])*log(pAB) +
			(2*n3x3[1,3]+n3x3[1,2]+n3x3[2,3])*log(pA-pAB) +
			(2*n3x3[3,1]+n3x3[2,1]+n3x3[3,2])*log(pB-pAB) +
			(2*n3x3[3,3]+n3x3[3,2]+n3x3[2,3])*log(1-pA-pB+pAB) + 
			n3x3[2,2]*log(pAB*(1-pA-pB+pAB) + (pA-pAB)*(pB-pAB))
	}
	
	
	solution <- optimize(
			loglik,
			lower=pmin + .Machine$double.eps,
			upper=pmax - .Machine$double.eps,
			maximum=TRUE
	)
	
	
	pAB <- solution$maximum
	
	estD <- pAB - pA*pB
	
	if (estD > 0)  
		estDp <- estD / Dmax
	else
		estDp <- estD / Dmin
	
	n <-  sum(n3x3)
	
	corr <- estD / sqrt( pA * pB * pa * pb )
	
	dchi <- (2*n*estD^2)/(pA * pa * pB* pb)
	dpval <- 1 - pchisq(dchi,1)
	
	retval <- list("r"=corr, "pAB"=pAB, "pA"=pA, "pB"=pB)
	#return r
	retval
}


main <- function()
{
	mleC <- readMLECounts()
	c11m <- mleC$c11m
	c12m <- mleC$c12m
	c13m <- mleC$c13m
	c21m <- mleC$c21m
	c22m <- mleC$c22m
	c23m <- mleC$c23m
	c31m <- mleC$c31m
	c32m <- mleC$c32m
	c33m <- mleC$c33m
	pA <- mleC$pA
	
	print(c11m)
	stop()
	
	m = nrow(c11m)
	n = ncol(c11m)
	
	mleR <- matrix(NA, m, n)
	
	
	for(i in 1:(m-1))
	{
		for(j in (i+1):m)
		{
			ret <- calculateRValue(pA[1,i], pA[1,j], 
					c11m[i,j], c12m[i,j], c13m[i,j],
					c21m[i,j], c22m[i,j], c23m[i,j],
					c31m[i,j], c32m[i,j], c33m[i,j])
			r <- ret$r
			mleR[i,j] <- mleR[j,i] <- r
		}
	}
	
	diag(mleR) <- NA
	
	print(mleR)
	
	stop()
	
	dumpR <- readDumpRValue()
	recoverRate <- singRecoverate(mleR, dumpR)
	
	cat("recover rate = ", recoverRate, "\n")
	
}

#test goes here
test <- function()
{
	ret <- calculateRValue(0.53825, 0.798, 223, 272, 71, 627, 394,0,413,0,0)
	correctValue = -0.4657458
	cat(ret$r, "=?", correctValue, "\n")
	
	ret <- calculateRValue(0.53825, 0.86175, 451,108,7,741,259,21,294,108,11)
	correctValue = 0.07892796
	cat(ret$r, "=?", correctValue, "\n")
	
	ret <- calculateRValue(0.53825, 0.56875, 177, 284,105,317,514,190,144,201,68)
	correctValue = -0.0268606
	cat(ret$r, "=?", correctValue, "\n")
	
	
#	if(ret$r == correctValue)
#	{
#		cat("test complete: \n")
#		print("PASS test\n")
#	}
}

debug(main)
main()







