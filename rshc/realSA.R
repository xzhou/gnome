# TODO: Add comment
# 
# Author: xzhou
###############################################################################


source("SHC.R")

#' using real R value for simulated annealing
#' 
#' @param var include some global information 
#' @param saConf is the configuration
realRSA <- function(var, saConf, ...)
{
	#debug(evaluate)
	print(var$targetGenoData)
	target <- list()
	target$RValues <- calculateRealR(var$targetGenoData)
	target$singleAlleleFreq <- calculateSingleAlleleFrequence(var$targetGenoData)
	
	cat("starting realRSA\n")
	finalPopList = list()
	for(ti in 1:saConf$totalIt)
	{
		#sink(file="sa.log")
		cat("start simulated annealing algorithm\n")
		T <- saConf$initT	#the init T
		t <- 1				#iteration counter
		i <- 1				#the iteration
		x <- generateRandomSample(var$nIndividuals, var$nSnps)
		cat(ncol(x), nrow(x), "\n")
		initSample <- list()
		initSample$genotype <- x
		initSample$RValues <- calculateRealR(x)
		initSample$singleAlleleFreq <- calculateSingleAlleleFrequence(x)
		currentQuality	<- evaluate(initSample, target)
		currentSample <- list()
		
		#print("init quality\n")
		#print(currentQuality)
		
		
		while(T >= saConf$Tmin && currentQuality$quality > 0)
		{
			for(i in 1:saConf$k)
			{
				newx <- singlePointMutate(x)
				sample <- list()
				sample$RValues <- calculateRealR(newx)
				sample$singleAlleleFreq <- calculateSingleAlleleFrequence(newx)
				sample$genotype <- newx
				newQuality <- evaluate(sample, target)
				
				if(newQuality$quality < currentQuality$quality)
				{
					x <- newx
					currentQuality <- newQuality
					currentSample <- sample
					cat(t, " - ", "qdiff = ", currentQuality$quality, "  signRecoverRate = ", currentQuality$recoverRate)
					cat("  ", "rdiff = ", currentQuality$normalizedRdiff, "  fdiff = ", currentQuality$normalizedFeqDiff, "\n")
					save(x, file = "currentPop")
					if(currentQuality$quality == 0)
					{
						#print(var$targetGenoData)
						#print(x)
						print(rbind(sortMatrixByRow(majorize(var$targetGenoData)), sortMatrixByRow(majorize(x))))
						save(x, file = "sampleGenoType")
						evaluate(sample = currentSample, target = target)	
						warning("target achieved\n")
						stop()
						break
					}		
				}
				#
				else
				{
					delta <-  newQuality$quality - currentQuality$quality
					p <- exp(-delta/T)
					randomV <- runif(1, 0, 1)
					if(randomV < p)
					{
						x <- newx
						currentQuality <- newQuality
						currentSample <- sample
						cat(t, " + ", "qdiff = ", currentQuality$quality, "  signRecoverRate = ", currentQuality$recoverRate)
						cat("  ", "rdiff = ", currentQuality$normalizedRdiff, "  fdiff = ", currentQuality$normalizedFeqDiff)
						cat("  ", "p = ", p, "\n")
						save(x, file = "currentPop")
					}
					else
					{
						#cat(t, "\tX\t Rej\n")
					}
				}
				t = t + 1
			}
			T <- saConf$beta*T #cool downc
			cat("cool down", T, "\n")
			
			if(currentQuality$quality < saConf$minDiff)
				break
		}
		fileName = paste("finalPop", ti, sep="")
		#finalPopList[i] = x
		save(x, file = fileName)
	}
	save(finalPopList, file = "finalPopList")
	cat("realRSA complete\n")
}


runRealRSA <- function()
{
	#-------------------------VARIABLES--------------------------------
	#var is the global configuration variable
	var <- list()
	saConf <- list()
	var$targetGenoData <- NULL
	var$max_it <- NULL
	var$nIndividuals <- NULL
	var$nSnps <- NULL
	
	#-------------------------START FROM HERE--------------------------
	#configuration
	var$max_it <- 1000000
	var$nIndividuals <- 2
	var$nSnps <- 10
	var$T <- 0.1	#for statistic hill climbing
	
	saConf$initT <- 0.01
	saConf$Tmin <- 1e-7	#minial temperature
	saConf$beta <- 0.9	#exponetial decreasing temperature
	saConf$k <- 10000		#number of iterations for each level of temperature
	saConf$totalIt = 10		#repeat the whole algorithm
	saConf$minDiff = 0.001
	
	cat("reading genodata from fasta file ...")
	var$targetGenoData <- readGenotypeFromFastaFile(nIndividuals = var$nIndividuals, nSnps = var$nSnps)
	targetGenoData <- var$targetGenoData
	save(targetGenoData, file = "targetGenoData")
	cat("complete \n")
	#print(var$targetGenoData)
	realRSA(var, saConf)
}
#debug(runRealRSA)
runRealRSA()
