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

testCalculateRealR()



