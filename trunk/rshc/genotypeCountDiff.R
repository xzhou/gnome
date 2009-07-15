# TODO: Add comment
# 
# Author: xzhou
###############################################################################

source("SHC.R")

calculateCounts <- function(genotype)
{
	genotype <- majorize(genotype)
	
	m <- nrow(genotype)
	n <- ncol(genotype)
	
	if(n%%2 != 0)
	{
		warning("not even snps data")
		stop()
	}
	
	nIndividuals <- m
	nSnps <- round(n/2)
	
	combinedGenotype <- matrix(0, 2*m, nSnps)
	
		#rearrange snps
	for(j in seq(1, n, 2))
	{
		left <- genotype[,j]
		right <- genotype[,j+1]
		
		colIndex = (j+1)/2
		
		for(i in seq(1, 2*m, 2))
		{
			rowIndex <- (i+1)/2
			combinedGenotype[i, colIndex] <- left[rowIndex]
			combinedGenotype[i+1, colIndex] <- right[rowIndex]
		}
	}
	
	c00 <- matrix(0.0, nSnps, nSnps)
	c00 <- c01 <- c02 <- c10 <- c11 <- c12 <- c20 <- c21 <- c22 <- c00
	
	for(i in seq(1, nSnps - 1))
	{
		for(j in seq(1, nSnps))
		{
			snp1 = combinedGenotype[,i]
			snp2 = combinedGenotype[,j]
			
			len <- length(snp1)
			
			for(k in seq(1, len, 2))
			{
				a1 = snp1[k]
				a2 = snp1[k+1]
				x1 = a1+a2
				b1 = snp2[k]
				b2 = snp2[k+1]
				x2 = b1+b2
				
				if(x1 == 2 && x2 == 2)
					c00[i,j] = c00[i,j] + 1
				if(x1 == 2 && x2 == 3)
					c01[i,j] = c01[i,j] + 1
				if(x1 == 2 && x2 == 4)
					c02[i,j] = c02[i,j] + 1
				if(x1 == 3 && x2 == 2)
					c10[i,j] = c10[i,j] + 1
				if(x1 == 3 && x2 == 3)
					c11[i,j] = c11[i,j] + 1
				if(x1 == 3 && x2 == 4)
					c12[i,j] = c12[i,j] + 1
				if(x1 == 4 && x2 == 2)
					c20[i,j] = c20[i,j] + 1
				if(x1 == 4 && x2 == 3)
					c21[i,j] = c21[i,j] + 1
				if(x1 == 4 && x2 == 4)
					c22[i,j] = c22[i,j] + 1
			}
		}
	}
	
	
	ret = list("c00" = c00,
			"c01" = c01,
			"c02" = c02,
			"c10" = c10,
			"c11" = c11,
			"c12" = c12,
			"c20" = c20,
			"c21" = c21,
			"c22" = c22
					)
	ret
}

testCountsDiff <- function()
{
	cat("reading genotype...")
	ceuGenotype = 
			readGenotypeFromFastaFile("/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/ceu.12encode")
	
	yriGenotype = 
			readGenotypeFromFastaFile("/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/yri.12encode")
	cat("complete\n")
	
	nCeuInd = nrow(ceuGenotype)
	nCeuSnp = ncol(ceuGenotype/2)
	
	nYriInd = nrow(yriGenotype)
	nYriSnp = ncol(yriGenotype/2)
	
	cat("nCeuInd = ", nCeuInd, "nCeuSnp = ", nCeuSnp, "nYriInd = ", nYriInd, "nYriSnp = ", nYriSnp, "\n")
	
	nInd = min(nCeuInd, nYriInd)
	nSnp = min(nCeuInd, nYriSnp)

	
		#select 
		ceuSample = ceuGenotype[1:nSnp, 1:(2*nSnp)]
		yriSample = yriGenotype[1:nSnp, 1:(2*nSnp)]
		
		cat("calculate ceuSample\n")
		#ceuCounts = calculateCounts(ceuSample)
		#save(ceuCounts, file = "ceuCounts")
		load("ceuCounts")
		
		cat("calculate yriSample\n")
		#yriCounts = calculateCounts(yriSample)
		#save(yriCounts, file = "yriCounts")
		load("yriCounts")
		
		upperIndex = upper.tri(ceuCounts$c00)
		
		#count diff
		d00 = (abs(ceuCounts$c00 - yriCounts$c00))[upperIndex]
		d01 = (abs(ceuCounts$c01 - yriCounts$c01))[upperIndex]
		d02 = (abs(ceuCounts$c02 - yriCounts$c02))[upperIndex]
		d10 = (abs(ceuCounts$c10 - yriCounts$c10))[upperIndex]
		d11 = (abs(ceuCounts$c11 - yriCounts$c11))[upperIndex]
		d12 = (abs(ceuCounts$c12 - yriCounts$c12))[upperIndex]
		d20 = (abs(ceuCounts$c20 - yriCounts$c20))[upperIndex]
		d21 = (abs(ceuCounts$c21 - yriCounts$c21))[upperIndex]
		d22 = (abs(ceuCounts$c22 - yriCounts$c22))[upperIndex]
		#print(d00)
		
		#percentage diff
		p00 = d00/(ceuCounts$c00[upperIndex])
		p01 = d01/(ceuCounts$c01[upperIndex])
		p02 = d02/(ceuCounts$c02[upperIndex])
		p10 = d10/(ceuCounts$c10[upperIndex])
		p11 = d11/(ceuCounts$c11[upperIndex])
		p12 = d12/(ceuCounts$c12[upperIndex])
		p20 = d20/(ceuCounts$c20[upperIndex])
		p21 = d21/(ceuCounts$c21[upperIndex])
		p22 = d22/(ceuCounts$c22[upperIndex])
		
		print(p00)

		
		diff = list("d00" = d00, "d01" = d01, "d02" = d02,
				"d10" = d10, "d10" = d11, "d02" = d12,
				"d20" = d20, "d21" = d21, "d22" = d22)
		
		#print diff
		meandiff = matrix(c(mean(d00), mean(d01), mean(d02), 
						mean(d10), mean(d11), mean(d12),
						mean(d20), mean(d21), mean(d22)), 3, 3, byrow = TRUE)
		
		maxdiff = matrix(c(max(d00), max(d01), max(d02), 
				max(d10), max(d11), max(d12),
				max(d20), max(d21), max(d22)), 3, 3, byrow = TRUE)

		mindiff = matrix(c(min(d00), min(d01), min(d02), 
				min(d10), min(d11), min(d12),
				min(d20), min(d21), min(d22)), 3, 3, byrow = TRUE)

		meanCounts = matrix(c(mean(ceuCounts$c00), mean(ceuCounts$c01), mean(ceuCounts$c02),
						mean(ceuCounts$c10), mean(ceuCounts$c11), mean(ceuCounts$c12),
						mean(ceuCounts$c20), mean(ceuCounts$c21), mean(ceuCounts$c22)), 3, 3, byrow = T)
		
		
		maxCounts =  matrix(c(max(ceuCounts$c00), max(ceuCounts$c01), max(ceuCounts$c02),
						max(ceuCounts$c10), max(ceuCounts$c11), max(ceuCounts$c12),
						max(ceuCounts$c20), max(ceuCounts$c21), max(ceuCounts$c22)), 3, 3, byrow = T)
		
		pmean = matrix(c(mean(p00), mean(p01), mean(p02), 
						mean(p10), mean(p11), mean(p12),
						mean(p20), mean(p21), mean(p22)), 3, 3)		
		
		
		#pmean = meandiff/meanCounts
		#pmax = maxdiff/maxCounts
		

		print(meandiff)
		print(maxdiff)
		
		print(meanCounts)
		print(maxCounts)
		
		print(pmean)
		#print(pmax)


}

testCountsDiff()

