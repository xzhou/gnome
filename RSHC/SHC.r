# TODO: Add comment
# 
# Author: xzhou
###############################################################################
source("rcalc.R")

shcMain <- function(targetGenotypeFileName = "", max_it = 1000, nIndividuals = 100, nSnps = 10)
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
		
		total <-  m*n*1.0
		correct <-  0.0
		for(i in 1:m)
			for(j in i:n)
			{
				if(targetGenotype[i,j] == sampleGenotype[i,j])
				{
					correct = correct + 1
				}
			}
		rate <- correct/total
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
				diff = targetRValues[i,j] - sampleRValues[i,j]
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

	
	stocasticHillClimb <- function(var, targetRValues, max_it, nIndividuals, nSnps)
	{
		x <- generateRandomSample(nIndividuals, nSnps)
		#print(x[1,])
		currentRValues = calculateRValues(x)
		#print(currentRValues)
		currentQuality <- evaluate(currentRValues, targetRValues)
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
				#sim <-  similarity(var.targetGenoData, newx)
				sim <- NULL
				cat(t, "\n RDiff=", currentQuality$diff,"\t signRecoverate", currentQuality$recoverRate, "\t similarity=", sim, "\n")
				#print(x)
				cat("\n\n")
			}
			else
			{
				cat(t, "\n")
			}
		}
	}
	
	#-------------------------VARIABLES--------------------------------
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
	
	cat("reading genodata from fasta file ...")
	var.targetGenoData <- readGenotypeFromFastaFile(nIndividuals = var.nIndividuals, nSnps = var.nSnps)
	cat("complete \n")
	cat("calculating real R value ... ... ...")
	var.targetRValue <- calculateRValues(var.targetGenoData)
	cat("complete ", nrow(var.targetRValue), "X", ncol(var.targetRValue), "\n")
	cat("\ntarget \n")
	#print(var.targetGenoData)
	#print(var.targetRValue[1,2:(var.nSnps)])
	stocasticHillClimb(var, var.targetRValue, var.max_it, var.nIndividuals, var.nSnps)
}

shcMain()
	