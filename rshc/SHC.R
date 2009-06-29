# TODO: Add comment
# 
# Author: xzhou
###############################################################################
source("rcalc.R")

#locate 
findGenotypeBlocks <- function(rValue, avgThreshold = 0.3, minBlockSize = 1)
{
	threshold = 0.6*avgThreshold
	rsValue <- rValue*rValue
	#print(rsValue)
	m <- nrow(rsValue)
	n <- ncol(rsValue)
	GenoBlocks <- matrix(0, 0, 4)
	i = 1
	while(i < m)
	{
		#cat(i, "\n")
		if(rsValue[i,i+1] >= threshold)
		{
			x = i			#row index
			y = i + 1		#col index
			#search horizontally
			while(rsValue[x,y] >= threshold)
			{
				y = y + 1
				if(y >n)
				{
					break
				}
			}
			yLim = y - 1
			xLim = y - 2
			
			#calculate the 
			total = 0.0
			count = 0
			for(x1 in x:xLim)
			{
				for(y1 in (x1+1):yLim)
				{
					count = count + 1
					total = total + rsValue[x1, y1]
				}
			}
			if((total/count) >= avgThreshold)
			{
				GenoBlocks = rbind(GenoBlocks, c(i, i+1, xLim, yLim))
			}
			i = xLim + 1
			cat("i = ", i, "xLim = ", xLim, "yLim = ", yLim, "\n")
		}
		else
		{
			i = i + 1
		}
	}
	#print(GenoBlocks)
	GenoBlocks
}

shcMain <- function(targetGenotypeFileName = "", max_it = 10000000, nIndividuals = 100, nSnps = 10)
{
	# read the file from ped file and convert it to standard genotype matrix
	readGenotypeFromFastaFile <- function(fileName = "../GenotypeLearnning/data/sim_4000seq/80SNP_CEU_sim_4000seq.12encode", nIndividuals = -1, nSnps = -1)
	{
		genoData <- read.table(fileName,
				header = FALSE)
		
		#DEBUG
		#print(genoData[1:2,])
		
		m <- nrow(genoData)
		n <- ncol(genoData)
		
		#DEBUG
		#cat(class(genoData), "nRow = ", m, "nCol = ", n, "\n")
		
		if(nIndividuals == -1)
			nIndividuals = m
		
		if(nSnps == -1)
			nSnps = n/2
		
		if(m < nIndividuals || n < 2*nSnps)
		{
			cat("not enough snps")
			stop()
		}
		
		#cut
		genoData <- genoData[1:nIndividuals, 1:(2*nSnps)]
		#as.character(genoData)
	}
	

	#generate a genotype matrix of nIndividuals X nSnps
	generateRandomSample <- function(nIndividual, nSnps)
	{
		x = rbinom(nIndividual*2*nSnps, 1, 0.5) + 1
		sampleGenoData <- matrix(x, nrow = nIndividual, ncol = 2*nSnps)
		as.data.frame(sampleGenoData)
	}
	
	#mutate a single point
	singlePointMutate <- function(genoData, ...)
	{
		m <- nrow(genoData)
		n <- ncol(genoData)
		
		i <- round(runif(1, 1, m))
		j <- round(runif(1, 1, n))
		
		if(genoData[i,j] == 1)
			genoData[i,j] = 2
		else if(genoData[i,j] == 2)
			genoData[i,j] = 1
		genoData	
	}
	
	
	#change 1 and 2 according to the 1 as major and 2 as minor
	majorize <- function(genoData, ...)
	{
		mGenoData <- NULL
		
		n <- ncol(genoData)
		m <- nrow(genoData)
		for(i in seq(1, n, by = 2))
		{
			snp1 = genoData[,i]
			snp2 = genoData[,i+1]
			
			combined = c(snp1, snp2)
			
			major.A = combined[which.max(combined)]
			major.a = combined[which.min(combined)]
			
			for (j in 1:m)
			{
				if(snp1[j] == major.A)
					snp1[j] = 1
				else
				{
					snp1[j] = 2
				}
				
				if(snp2[j] == major.A)
				{
					snp2[j] = 1
				}
				else
				{
					snp2[j] = 2
				}
			}
			mGenoData <- cbind(mGenoData,  snp1, snp2)	
		}
		names(mGenoData) = NULL
		mGenoData
	}
	
	#calculate the R values
	calculateRValues <- function(genoData)
	{
		formatedGenotype <- NULL
		n <-  ncol(genoData)
		
		#DEBUG
		#cat("genoData", nrow(genoData), ncol(genoData), "\n")
		
		for(i in seq(1, n, by = 2))
		{
			formatedGenotype = cbind(formatedGenotype, cbind(paste(genoData[,i], genoData[,i+1], sep="/")))
		}
		
		#debug
		#print(formatedGenotype[,77])
		
		n <-  ncol(formatedGenotype)
		m <-  nrow(formatedGenotype)
		
		#cat("formated", m,n)
		
		plotGenoData <- data.frame(formatedGenotype)
		
		#DEBUG
		#cat("plotGenoData", nrow(plotGenoData), ncol(plotGenoData), "\n")
		
		for(i in seq(1, n))
		{
			plotGenoData[,i] <- genotype(formatedGenotype[,i])
		}
		
		#print(plotGenoData[,77])
		retValue <- calcAllR(plotGenoData)
		
		#DEBUG
		#stop("Debug check @calculateRValue")
		
		rValue <- retValue$r
	}
	
	sortMatrixByRow <- function(aGenotype, ...)
	{
		m = nrow(aGenotype)
		n = ncol(aGenotype)
		
		more <- function(a, b)
		{
			len = length(a)
			for(i in 1:len)
			{
				if(a[i] > b[i])
				{
					return(TRUE)
				}
				else if(a[i] < b[i])
				{
					return(FALSE)
				}
			}
			#the same
			return(FALSE)
		}
		
		#bubble sort
		swapped = FALSE
		for(i in 1:(m-1))
		{
			for (j in 1:(m-i))
			{				
				#print(aGenotype[j,])
				#print(aGenotype[j+1, ])
				#cat("j = ", j, "\n")
				if(more(aGenotype[j,], aGenotype[j+1,]))
				{
					temp = aGenotype[j, ]
					aGenotype[j,] = aGenotype[j+1, ]
					aGenotype[j+1, ] = temp
					swapped = TRUE
				}
			}
			if(!swapped)
				break
		}
		aGenotype
	}
	
	similarity <- function(targetGenotype, sampleGenotype)
	{
		#majorize both of the genotype
		targetGenotype <- majorize(targetGenotype)
		sampleGenotype <- majorize(sampleGenotype)
		
		m <- mt <-  nrow(targetGenotype)
		n <- nt <-  ncol(sampleGenotype)
		
		ms <-  nrow(sampleGenotype)
		ns <-  ncol(sampleGenotype)
		
		if( mt != ms || ns != ns)
		{
			warning("incompatible geno type")
			stop()
		}
		#sort the genotype to compare
		targetGenotype = sortMatrixByRow(targetGenotype)
		sampleGenotype = sortMatrixByRow(sampleGenotype)
		
		totalElement = m*n
		correct = 0.0
		for(i in 1:m)
		{
			for(j in 1:n)
			{
				if(targetGenotype[i,j] == sampleGenotype[i,j])
				{
					correct = correct + 1
				}
			}
		}
		
		#return the rate
		correct/totalElement
	}
	
	
	#evaluate the population
	evaluate <- function(sampleRValues, targetRValues, ...) 
	{
		#sampleRValues <- calculateRValues(genoData)
		totalDifference = 0.0
		totalSigns = 0.0
		correctSings = 0.0
		m = nrow(sampleRValues)
		n = ncol(sampleRValues)
		
		if( m != n)
		{
			warning("different matrix")
			stop()
		}
		#print(sampleRValues)
		#print(targetRValues)
		for(i in 1:m)
		{
			for(j in 1:n)
			{
				#calculate the R square difference as the 
				diff = targetRValues[i,j]*targetRValues[i,j] - sampleRValues[i,j]*sampleRValues[i,j]
				totalDifference = totalDifference + abs(diff)
				totalSigns = totalSigns + 1
				if((targetRValues[i,j]*sampleRValues[i,j]) >= 0)
				{
					correctSings = correctSings + 1
				}
			}
		}
		
		ret <- list("diff"=totalDifference, "recoverRate" = correctSings/totalSigns)
	}

	#simulated annealing algorithm 
	sa <- function(var, saConf, ...)
	{
		for(ti in 1:saConf.totalIt)
		{
			#sink(file="sa.log")
			cat("start simulated annealing algorithm\n")
			T <- saConf.initT	#the init T
			t <- 1				#iteration counter
			i <- 1				#the iteration
			x <- generateRandomSample(var.nIndividuals, var.nSnps)
			cat(ncol(x), nrow(x), "\n")
			currentRValues <- calculateRValues(x)
			currentQuality	<- evaluate(currentRValues, var.targetRValue)
			
			while(T >= saConf.Tmin)
			{
				for(i in 1:saConf.k)
				{
					newx <- singlePointMutate(x)
					newRValues <- calculateRValues(newx)
					newQuality <- evaluate(newRValues, var.targetRValue)
					
					if(newQuality$diff < currentQuality$diff)
					{
						x <- newx
						currentRValue <- newRValues
						currentQuality <- newQuality
						cat(t, "\t-\t", "diff = ", currentQuality$diff, "\tsignRecoverRate = ", currentQuality$recoverRate, "\n")
						save(x, file = "currentPop")
					}
					else
					{
						delta <-  newQuality$diff - currentQuality$diff
						p <- exp(-delta/T)
						randomV <- runif(1, 0, 1)
						if(randomV < p)
						{
							x <- newx
							currentRValue <- newRValues
							currentQuality <- newQuality
							cat(t, "\t+\t", "diff = ", currentQuality$diff, "\tsignRecoverRate = ", currentQuality$recoverRate, "\t", p, "\n")
							save(x, file = "currentPop")
						}
						else
						{
							cat(t, "\tX\t Rej\n")
						}
					}
					t = t + 1
				}
				T <- saConf.beta*T #cool downc
				cat("cool down", T, "\n")
			}
			fileName = paste("finalPop", ti, sep="")
			save(x, file = fileName)
		}
	}
	
	#stacastic algorithm
	stocasticHillClim <- function(var,...)
	{
		cat("start stocastic hill climbing with max_it = ", var.max_it, "T = ", var.T, "\n")
		x <- generateRandomSample(var.nIndividuals, var.nSnps)
		currentRValues <- calculateRValues(x)	#get the R values
		
		currentQuality <- evaluate(currentRValues, var.targetRValue)
		
		#print(currentQuality)
		
		t <- 0
		
		while(t < var.max_it && currentQuality$diff != 0)
		{
			t <-  t + 1
			newx <- singlePointMutate(x)
			newRValues <- calculateRValues(newx)
			newQuality <- evaluate(newRValues, var.targetRValue)
			
			diff <-  newQuality$diff - currentQuality$diff
			p <- 1/(1+exp(diff/var.T))
			aRandomNumber = runif(1, 0, 1)			
			if(aRandomNumber  < p){
				x <- newx
				currentQuality <-  newQuality
				currentRValues <- newRValues
				cat(t, "RDiff=", currentQuality$diff,"\t signRecoverate", currentQuality$recoverRate, "\n")
			}
			else
			{
				cat(t, "\n")
				#cat(t, "p =", p)
			}
		}
		
	}
	
	hillClimb <- function(var, targetRValues, max_it, nIndividuals, nSnps)
	{
		x <- generateRandomSample(nIndividuals, nSnps)
		#print(x[1,])
		currentRValues = calculateRValues(x)
		#print(currentRValues)
		currentQuality <- evaluate(currentRValues, targetRValues)
		cat("0 RDiff =", currentQuality$diff, "\t signRecoverate =", currentQuality$recoverRate,"\n")
		t <- 0
		while (t < max_it && currentQuality$diff != 0 && currentQuality$recoverRate != 1)
		{
			t <- t + 1
			newx <- singlePointMutate(x)
			newRValues <- calculateRValues(newx)
			newQuality <- evaluate(newRValues, targetRValues)
			if (newQuality$diff < currentQuality$diff)
			{
				x <- newx
				currentQuality <-  newQuality
				currentRValues <- newRValues
				sim <-  similarity(var.targetGenoData, newx)
				cat(t, " RDiff=", currentQuality$diff,"\t signRecoverate", currentQuality$recoverRate, "\t similarity=", sim, "\n")
				#print(x)
			}
			else
			{
				cat(t, "\n")
			}
		}
	}
	
	#-------------------------VARIABLES--------------------------------
	#var is the configuration variable
	var.targetGenoData <- NULL
	var.max_it <- NULL
	var.nIndividuals <- NULL
	var.nSnps <- NULL
	var.currentGenoMatrxi <- NULL	
	
	#-------------------------START FROM HERE--------------------------
	#init variables	
	var.max_it <- max_it
	var.nIndividuals <- nIndividuals
	var.nSnps <- nSnps
	
	var.T <- 0.1	#for statistic hill climbing
	
	saConf.initT <- 0.1
	saConf.Tmin <- 0.001	#minial temperature
	saConf.beta <- 0.8	#exponetial decreasing temperature
	saConf.k <- 100		#number of iterations for each level of temperature
	saConf.totalIt = 10		#repeat sa algorithm for times to see if you have 
							#multiple local optimal value
	
	
	cat("reading genodata from fasta file ...")
	var.targetGenoData <- readGenotypeFromFastaFile(nIndividuals = var.nIndividuals, nSnps = var.nSnps)
	cat("complete \n")
	cat("calculating real R value ... ... ...")
	#var.targetRValue <- calculateRValues(var.targetGenoData)
	load("rValue")
	cat("complete ", nrow(var.targetRValue), "X", ncol(var.targetRValue), "\n")
    
	#stocasticHillClim(var, var.targetRValue, var.max_it, var.nIndividuals, var.nSnps)
	#debug(sa)
	sa(var, saConf)
}
#run the function as default config
#debug(findGenotypeBlocks)
shcMain()