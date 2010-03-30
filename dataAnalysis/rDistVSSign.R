# TODO: Analysis the relationship between sing and r valude distribtuion
# 
# Author: xzhou
###############################################################################

source("SHC.R")

rDistVSSign <- function(totalGenotype, it = 100)
{
	m <- nrow(totalGenotype)
	n <- ncol(totalGenotype)
	
	#number individual sampled
	sampleIndividualSize = round(m/10)
	
	#total r values
	totalRVs <- calculateRValues(totalGenotype)
	
	nSnps <- n/2
	
	result <- NULL
	
	
	#iterate through each pair
	for(i in 1:(nSnps-1))
	{
		for(j in (i+1):nSnps)
		{
			
			pairGenotype <- totalGenotype[, c(2*i-1, 2*i, 2*j-1, 2*j)]
			
			totalV <- totalRVs[i,j]
			
			aResult <- c(i,j,totalV)
			
			for(k in 1:it)
			{
				#randomly select sampleIndividualSize people
				pool = 1:m
				sampleIndex <- sort(sample(pool, sampleIndividualSize))
				sampledGenotype <- pairGenotype[sampleIndex, ]
				rValue <- calculateRValues(sampledGenotype)
				aResult <- c(aResult, rValue[1,2])
			}
			#print(aResult)
			
			hist(aResult[-1:-3], xlab = paste(aResult[1]), ylab = paste(aResult[2]), breaks = 20)
			abline(v = aResult[3], col="blue", lty = 1, lwd = 3)
			result <- rbind(result, aResult)
		}
	}
	
	save(file = "rdist.RData", result)
	result
}

rdist <- function(file = "../data/sim_4000seq/80SNP_CEU_sim_4000seq.12encode", popSize = 2000, nSnps = 30 )
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
	
	pdf("dist.pdf")
	result <- rDistVSSign(targetGenotype)
	dev.off()
}

dist2Sign <- function(result = NULL, loadFromFile = F)
{
	
	#DEBUG
	#loadFromFile = T
	
	if(result == NULL || loadFromFile == T)
	{
		load("rdist.RData")
	}
	
	nSnps <-  max(result[,2])
	
	signMatrix = matrix(0, nSnps, nSnps)
	
	nRec <- nrow(result)
	
	for(k in 1:nRec)
	{
		aResult <- result[k, ]
		i <- aResult[1]
		j <- aResult[2]
		
		realR <- aResult[3]
		
		allAvgR <- aResult[-1:-3]
		
		nPositive <- length(allAvgR[allAvgR>0])
		nNegative <- length(allAvgR[allAvgR<0])
		
		if(nPositive > nNegative)
		{
			signMatrix[i,j] <- signMatrix[j,i] <-  1
		}
		else if(nPositive < nNegative)
		{
			signMatrix[i,j] <- signMatrix[j,i] <-  -1
		}
		else
		{
			signMatrix[i,j] <-  0
		}
	}
	signMatrix
}

distPlot <- function(result = NULL, loadFromFile = T)
{
	
	if(loadFromFile == TRUE || result == NULL)
	{
		load("rdist.RData")
	}
	
	m <- nrow(result)
	n <- ncol(result)
	pdf("adjust_r_dist_2000_10_100__0.05.pdf")
	for(i in 1:m)
	{
		aResult = result[i,]
		#hist(aResult[-1:-3], xlab = paste(aResult[1]), ylab = paste(aResult[2]), breaks = 20)
		if(abs(aResult[3]) >= 0.05)
		{
			hist(aResult[-1:-3], xlab = paste(aResult[1]), ylab = paste(aResult[2]), breaks = seq(-1, 1, by = 0.05), xlim = c(-1, 1), main = paste(aResult[1], "-", aResult[2], "r = ", aResult[3]))
			abline(v = aResult[3], col="blue", lty = 1, lwd = 3)
			abline(v = 0.05, col="red", lty = 3, lwd = 1)
			abline(v = -0.05, col="red", lty = 3, lwd = 1)
		}
	}
	dev.off()
}

#rdist()
sign <- dist2Sign()
print(sign)

