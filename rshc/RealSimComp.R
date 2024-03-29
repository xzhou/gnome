# TODO: Add comment
# 
# Author: xzhou
###############################################################################

source("SHC.R")

t <- function()
{
	t <- readGenotypeFromFastaFile()
	
	er <- calculateRValues(t)
	save(er, file = "77_2000estR")
	
	rr <- calculateRealR(t)
	save(rr, file = "77_2000realR")
	
	sr <- singRecoverate(er, rr)
	
	print(sr)
}

t()