# TODO: analysis how individuals we need to calculate get the correct sign
# 
# Author: xzhou
###############################################################################

source("SHC.R")

popSizeVSSignRate <- function(totalGenotype, it = 5)
{
	m <- nrow(totalGenotype)
	n <- ncol(totalGenotype)
	
	rValue <- calculateRValues(totalGenotype)
	
	if(m < 50)
	{
		warning("sample too small, we need at least 100 individual\n")
		stop()
	}
	
	#the result of k and sign recover rate
	result <- NULL
	
	for(k in seq(10, m, by = 50))
	{
		
		totalSignRecoverRate = 0.0
		for(i in 1:it)
		{
			#randomly select k individual from totalGenotype and calculate the R value
			pool <- seq(1:m)
			sampleRows <- sort(sample(pool, k))
			
			sampleGenotype <- totalGenotype[sampleRows, ]
			
			sampleR <- calculateRValues(sampleGenotype)
			
			sampleRecoverRate <- singRecoverate(sampleR, rValue)
			
			totalSignRecoverRate = totalSignRecoverRate + sampleRecoverRate
		}
		
		avgRecoverate = totalSignRecoverRate/it
		
		cat(k, avgRecoverate, "\n")
		
		result = rbind(result, c(k, avgRecoverate))
	}
	
	
	#draw relationship
	save(file = "result.RData", result)
	plot(result[,1], result[,2])
	dev.print(device = pdf, file = "popSizeVsSignRate.pdf")
}

Entry <- function(file = "../data/sim_4000seq/80SNP_CEU_sim_4000seq.12encode", popSize = -1, nSnps = -1)
{
	targetGenotype <- read.table(file)
	
	m <- nrow(targetGenotype)
	n <- ncol(targetGenotype)
	
	if(popSize == -1)
	{
		popSize = m
	}
	if(nSnps == -1)
	{
		nSnps = n/2
	}
	
	targetGenotype = targetGenotype[1:popSize, 1:(2*nSnps)]
	
	popSizeVSSignRate(targetGenotype)
}

Entry()



