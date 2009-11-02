# TODO: Add comment
# 
# Author: xzhou
###############################################################################

source("SHC.R")

test <- function()
{
	genotype = matrix(0, 4, 4)
	genotype[1,] = c(NA, 1, 1, 2)
	genotype[2,] = c(1, 1, 2, 2)
	genotype[3,] = c(1, 1, 2, 1)
	genotype[4,] = c(1, 1, 1, 2)
	
	genotype <- as.data.frame(genotype)
	
	singleAlleleFreq <- calculateSingleAlleleFrequence(genotype)
	
	print(singleAlleleFreq)
}

test()