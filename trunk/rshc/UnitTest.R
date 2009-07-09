# TODO: This is a unit test module for all the function of SHC.R
# 
# Author: xzhou
###############################################################################

source("SHC.R")

testCalculateRealR <- function()
{
	targetGenotype <- readGenotypeFromFastaFile()
	
	rValue <- calculateRealR(targetGenotype)
	
	write.table(rValue, file = "realRValue")
	#compare result
	rValue
}

testFreq <- function()
{
	genotype <- matrix(round(runif(20, 1, 2)), 2, 10)
}



equalTest <- function()
{
	x1 		= matrix(c(2,1,1,2,1,1,2,1,1,2,1,1,1,1,2,2,2,1,2,1,1,1,2,2,1,1,2,2,1,1,1,2,1,2,1,1,1,1,1,1), nrow = 2, byrow = T)
	x2 		= matrix(c(1,1,2,2,1,1,2,1,2,1,2,2,1,1,1,2,2,2,2,2,2,1,1,2,1,1,1,1,1,1,2,1,1,2,2,1,1,2,1,2), nrow=2, byrow = T)
	target 	= matrix(c(2,1,1,2,1,1,1,1,2,2,2,1,1,2,2,2,1,2,2,1,1,1,2,2,1,1,1,2,2,1,2,2,1,1,1,2,1,2,2,1), nrow = 2, byrow = T)
}

#testCalculateRealR()



