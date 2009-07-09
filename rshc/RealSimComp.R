# TODO: Add comment
# 
# Author: xzhou
###############################################################################

source("SHC.R")

t <- function()
{
	t <- readGenotypeFromFastaFile()
	
	er <- calculateRValues(t)
	rr <- calculateRealR(t)
	
	singRecoverate(er, rr)
}

t()