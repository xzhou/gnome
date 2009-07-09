# TODO: Add comment
# 
# Author: xzhou
###############################################################################

source("SHC.R")

t <- function()
{
	t <- readGenotypeFromFastaFile()
	
	er <- calculateRValues(t)
	save(er, "77_2000estR")
	
	rr <- calculateRealR(t)
	save(rr, "77_2000realR")
	
	singRecoverate(er, rr)
}

t()