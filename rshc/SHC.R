# TODO: Add comment
# 
# Author: xzhou
###############################################################################
source("rcalc.R")

#change 1 and 2 according to the 1 as major and 2 as minor
majorize <- function(genoData, verbose = F, ...)
{
	mGenoData <- NULL
	
	n <- ncol(genoData)
	m <- nrow(genoData)
	for(i in seq(1, n, by = 2))
	{
		snp1 = genoData[,i]
		snp2 = genoData[,i+1]
		
		combined = c(snp1, snp2)
		c1 = length(combined[combined == 1])
		c2 = length(combined[combined == 2])
		
		if(c1 >= c2 )
		{
			major.A = 1
			major.a = 2
		}
		else
		{
			major.A = 2
			major.a = 1
		}	
#		
#		#major.A = combined[which.max(combined)]
#		#major.a = combined[which.min(combined)]
#		fct <- as.factor(combined)
#		smr <- summary(fct[!is.na(fct)])
#		
#		if(length(smr) == 1)
#		{
#			major.A = 1
#			minor.A = 2
#		}
#		else if(length(smr) == 2)
#		{
#			major.A = as.numeric(names(which.max(smr)))
#			minor.A = as.numeric(names(which.min(smr)))
#		}
#		else
#		{
#			warning("each allele should have 2 types of nucleotide")
#			stop()
#		}
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

#calculate real R value using haplotype to see if any imporovements
calculateRealR <- function(genotype, verbose = F)
{
	#calculare major
	genotype <- majorize(genotype)
	
	m <- nrow(genotype)
	n <- ncol(genotype)
	
	nIndividuals <- m
	nSnps <- round(n/2)
	
	combinedGenotyp <- matrix(0, 2*m, nSnps)
	
	#rearrange snps
	for(j in seq(1, n, 2))
	{
		left <- genotype[,j]
		right <- genotype[,j+1]
		
		colIndex = (j+1)/2
		
		for(i in seq(1, 2*m, 2))
		{
			rowIndex <- (i+1)/2
			combinedGenotyp[i, colIndex] <- left[rowIndex]
			combinedGenotyp[i+1, colIndex] <- right[rowIndex]
		}
	}
	
	#print(combinedGenotyp)
	m <- nrow(combinedGenotyp)
	n <- ncol(combinedGenotyp)
	
	r <- matrix(-2, n, n)
	
	diag(r) <- 0
	
	p00 <- r
	
	for(i in seq(1, n-1))
	{
		for(j in (i+1):n)
		{
			c00 <- c01 <- c10 <- c11 <- 0
			
			snp1 = combinedGenotyp[,i]
			snp2 = combinedGenotyp[,j]
			
			len <- length(snp1)
			
			#modify the code for
			for(k in seq(1, len))
			{
				if(snp1[k] == 1 && snp2[k] == 1)
					c00 = c00 + 1
				else if(snp1[k] == 1 && snp2[k] == 2)
					c01 = c01 + 1
				else if(snp1[k] == 2 && snp2[k] == 1)
					c10 = c10 + 1
				else if(snp1[k] == 2 && snp2[k] == 2)
					c11 = c11 + 1
				else
				{
					warning("error reading genotype")
					stop()
				}
			}
			
			pAB = c00/(c00+c01+c10+c11)
			p00[i,j] <- p00[j,i] <-  pAB
						
			c0x <- c00 + c01
			c1x <- c10 + c11
			cx0 <- c00 + c10
			cx1 <- c01 + c11
			
			D = c00*c11 - c01*c10 + 0.0
			L = c0x*c1x*cx0*cx1
			
			if(L == 0)
			{
				r[j,i] <- r[i,j] <- 0
			}
			else
			{
				r[j,i] <- r[i,j] <- D/sqrt(L)
				#cat(i,j, r[i,j], "\n")
			}
		}
	}
	
	if(verbose == T)
	{
		print(r)
	}
	
	diag(r) <- NA
	diag(p00) <- NA
	retL <- list("r" = r, "p00"=p00)
	retL
}

diffR <- function(estR, realR)
{
	fileName <- "/home/xzhou/research_linux/gnome/workspace/GenotypeLearnning/data/sim_4000seq/80SNP_CEU_sim_4000seq.rValue"
	realR <- read.table(filename)
}


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

#it is difficult to compare the distance between to genotype since the 
#individual is not aligned, we use a greedy algorithm to determine the distance

maxWeightSimilarity <- function(g1, g2, ...)
{
	m1 = nrow(g1)
	n1 = ncol(g1)
	m2 = nrow(g2)
	n2 = ncol(g2)
	
	if(m1 != m2 || n1 != n2)
	{
		warning("inconsistent matrix")
		stop()
	}
	
	similarityMatrix = matrix(0.0, m1, m2)
	#calculate the distance of each row pairwise
	for(i in 1:m1)
	{
		for(j in 1:m2)
		{
			r1 = g1[i,]
			r2 = g2[j,]
			diff = r1 - r2
			sim = length(which(diff == 0))
			similarityMatrix[i,j] = sim
		}
	}
	
	m = similarityMatrix
	
	#greedy algorithm
	totalSim = 0.0
	for(i in 1:m1)
	{
		x = which.max(m)

		colIndex = ceiling(x/m1)
		rowIndex = x - (m1*(colIndex-1))
		
		#cat("<", colIndex, rowIndex, "> max = ", m[x], "\n")
		
		totalSim = m[x] + totalSim
		
		m[,colIndex] = -1
		m[rowIndex, ] = -1

	}
	totalSim
}

checkSimlarity <- function()
{
	crossCompare <- matrix(0.0, 10, 10)
	popList <- list()
	for(i in 1:9)
	{
		fileName = paste("finalPop", i, sep = "")
		load(fileName)
		y = x
		for(j in (i+1):10)
		{
			fileName = paste("finalPop", j, sep = "")
			load(fileName)
			z = x

			sim = maxWeightSimilarity(y, z)
			
			crossCompare[i,j] = sim
			crossCompare[j,i] = sim
			cat("<", i, j, "> = ", sim, "\n")
		}
	}
	
	print(crossCompare)
	
}

convert <- function()
{
	for(i in 1:10)
	{
		fileName = paste("finalPop", i, sep="")
		load(fileName)
		write.table(x, file = paste(fileName, ".txt", sep=""), col.names=FALSE, row.names=FALSE)
	}
}

formateGenotype <- function(genoData)
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
	
	plotGenoData
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


#calculate single allele frequency
calculateSingleAlleleFrequence <- function(genotype, ...)
{
	m <- nrow(genotype)
	n <- ncol(genotype)
	
	if(n%%2 != 0)
	{
		warning("genotype should have even number of alleles")
		stop()
	}
	
	nSnps <- round(n/2)
	
	#first row = major counts
	#second row = minor counts
	singleAlleleFrequency = matrix(-1, nrow = 2, ncol = nSnps)
	
	for(i in 1:nSnps)
	{
		snp <- c(genotype[,2*i-1], genotype[,2*i])
		fct <- as.factor(snp)
		smr <- summary(fct[!is.na(fct)])
		
		if(length(smr) == 1)
		{
			majorCount <- max(smr)
			minorCount <- 0
		}
		else if(length(smr) == 2)
		{
			majorCount <- max(smr)
			minorCount <- min(smr)
		}
		else
		{
			warning("each allele should have only ")
			stop()
		}
		
		singleAlleleFrequency[1, i] <- majorCount
		singleAlleleFrequency[2, i] <- minorCount
	}
	
	singleAlleleFrequency
}


singRecoverate <- function(sampleRValues, targetRValues, verbose = F)
{
	m = nrow(sampleRValues)
	n = ncol(targetRValues)
	
	totalSigns = 0.0
	correctSigns = 0.0
	
	for(i in 1:(m-1))
	{
		for(j in (i+1):n)
		{
			totalSigns = totalSigns + 1
			if(targetRValues[i,j] != 0 && sampleRValues[i,j] != 0)
			{
				if(targetRValues[i,j] * sampleRValues[i,j] > 0)
				{
					correctSigns = correctSigns + 1
				}
			}
			else if(targetRValues[i,j] == 0 && sampleRValues[i,j] == 0)
			{
				correctSigns = correctSigns + 1
			}
		}
	}
	signRecoverRate <- correctSigns/totalSigns
}


#' Evaluate a sample
#' 
#' @param sample the sample genotype with sample.RValues, sample.singleAlleleFreq
#' @param target the target
evaluate <- function(sample = NULL, target = NULL, rWeight = 0.7, freqWeight = 0.3, verbose = F) 
{
	sampleRValues <- sample$RValues
	sampleSingleAlleleFreq <- sample$singleAlleleFreq
	
	targetRValues <- target$RValues
	targetSingleAlleleFreq <- target$singleAlleleFreq
	
	if(verbose)
	{
		cat("sample r values \n")
		print(sampleRValues)
		
		cat("real r values \n")
		print(targetRValues)
	}
	
	totalSigns = 0.0
	correctSigns = 0.0
	m = nrow(sampleRValues)
	n = ncol(sampleRValues)
	
	if( m != n)
	{
		warning("different matrix")
		stop()
	}
	#print(sampleRValues)
	#print(targetRValues)
	
	#upper tri angle
	
	for(i in 1:(m-1))
	{
		for(j in (i+1):n)
		{
			totalSigns = totalSigns + 1
			if(targetRValues[i,j] != 0 && sampleRValues[i,j] != 0)
			{
				if(targetRValues[i,j] * sampleRValues[i,j] > 0)
				{
					correctSigns = correctSigns + 1
				}
			}
			else if(targetRValues[i,j] == 0 && sampleRValues[i,j] == 0)
			{
				correctSigns = correctSigns + 1
			}
		}
	}
	signRecoverRate <- correctSigns/totalSigns
	
	#calc R diff
	rdiff = sum(abs(targetRValues*targetRValues - sampleRValues*sampleRValues))
	normalizedRDiff <- rdiff/sum(targetRValues*targetRValues)
	
	freqDiff = sum(abs(sampleSingleAlleleFreq - targetSingleAlleleFreq))
	normalizedFreqDiff = freqDiff*1.0/sum(targetSingleAlleleFreq)
	
	quality = rWeight*normalizedRDiff + freqWeight * normalizedFreqDiff
	
	ret <- list("quality" = quality, "normalizedRdiff"=normalizedRDiff, "normalizedFeqDiff" = normalizedFreqDiff, "recoverRate" = signRecoverRate)
}


#comfirm the result
checkFinalPop <- function()
{
	load("rValue")				#load var.targetRValue
	load("targetGenoData")		#load var.targetGenoData
	
	for(i in 1:10)
	{
		fileName <- paste("finalPop", i, sep = "")
		load(fileName)			#load x
		sampleR <- calculateRValues(x)
		result <- evaluate(sampleR, var.targetRValue)
		cat(i, "totalDiff = ", result$diff, "signRecoverRate", result$recoverRate, "\n")
	}
}

# read the file from ped file and convert it to standard genotype matrix
readGenotypeFromFastaFile <- function(fileName = "../data/sim_4000seq/80SNP_CEU_sim_4000seq.12encode", nIndividuals = -1, nSnps = -1)
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
	genoData
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

#TODO change code to accomodate the new evluate function
#simulated annealing algorithm 
sa <- function(var, saConf, ...)
{
	finalPopList = list()
	for(ti in 1:saConf.totalIt)
	{
		#sink(file="sa.log")
		cat("start simulated annealing algorithm\n")
		T <- saConf.initT	#the init T
		t <- 1				#iteration counter
		i <- 1				#the iteration
		x <- generateRandomSample(var.nIndividuals, var.nSnps)
		cat(ncol(x), nrow(x), "\n")
		#currentRValues <- calculateRValues(x)
		currentQuality	<- evaluate(currentRValues, var.targetRValue)
		
		while(T >= saConf.Tmin && currentQuality$diff > 0)
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
					save(x, file = "currentPop0diff")
					if(currentQuality$diff == 0)
						break
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
			
			if(currentQuality$diff < saConf.minDiff)
				break
		}
		fileName = paste("finalPop", ti, sep="")
		#finalPopList[i] = x
		save(x, file = fileName)
	}
	save(finalPopList, file = "finalPopList")
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

shcMain <- function(targetGenotypeFileName = "../data/sim_4000seq/80SNP_CEU_sim_4000seq.12encode", max_it = 10000000, nIndividuals = 100, nSnps = 10)
{	
	#-------------------------VARIABLES--------------------------------
	#var is the global configuration variable
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
	
	saConf.initT <- 0.01
	saConf.Tmin <- 1e-7	#minial temperature
	saConf.beta <- 0.8	#exponetial decreasing temperature
	saConf.k <- 100		#number of iterations for each level of temperature
	saConf.totalIt = 10		#repeat sa algorithm for times to see if you have 
							#multiple local optimal value
	saConf.minDiff = 0.001
	
	cat("reading genodata from fasta file ...")
	var.targetGenoData <- readGenotypeFromFastaFile(nIndividuals = var.nIndividuals, nSnps = var.nSnps)

	save(var.targetGenoData, file = "targetGenoData")
	print(var.targetGenoData)
	cat("complete \n")
	#realRSA(var, saConf)
}
