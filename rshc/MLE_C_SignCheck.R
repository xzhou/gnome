# Check if the sign of MLE recovered C has a high sign recover rate
# 
# Author: xzhou
###############################################################################

# read the dump R value
readDumpRValue <- function()
{
	
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

debug(readMLECounts)
readMLECounts(T)


