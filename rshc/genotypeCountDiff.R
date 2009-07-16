# Using reference group genotype count to reduce the searching space of 
# genotype counts
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
		for(j in seq(i+1, nSnps))
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
	
	for(i in 1:(nSnps-1))
	{
		for(j in (i+1):nSnps)
		{
			c00[j,i] = c00[i,j]
			c01[j,i] = c01[i,j]
			c02[j,i] = c02[i,j]
			c10[j,i] = c10[i,j]
			c11[j,i] = c11[i,j]
			c12[j,i] = c12[i,j]
			c20[j,i] = c20[i,j]
			c21[j,i] = c21[i,j]
			c22[j,i] = c22[i,j]
		}
	}
	
	diag(c00) <- diag(c01) <- diag(c02) <- diag(c10) <- diag(c11) <- diag(c12) <- diag(c20) <- diag(c21) <- diag(c22) <- NA
	
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

testCountsDiff <- function(loadSave = TRUE)
{
	cat("reading genotype...")
	ceuGenotype = 
			readGenotypeFromFastaFile("/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/ceu.12encode")
	
	yriGenotype = 
			readGenotypeFromFastaFile("/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/yri.12encode")
	cat("complete\n")
	
	
	ceur <- calculateRealR(ceuGenotype)
	yrir <- calculateRealR(yriGenotype)
	
	write.table(ceur, file = "ceur.txt")
	write.table(yrir, file = "yrir.txt")
	
	yris = yrir*yrir
	ceus = ceur*ceur
	
	yrisL = yris[yris>0.01]
	ceusL = ceus[ceus>0.01]
	
	hist(diff)
	dev.print(device=pdf, "rdiff001.pdf")
	
	nCeuInd = nrow(ceuGenotype)
	nCeuSnp = ncol(ceuGenotype)/2
	
	nYriInd = nrow(yriGenotype)
	nYriSnp = ncol(yriGenotype)/2
	
	cat("nCeuInd = ", nCeuInd, "nCeuSnp = ", nCeuSnp, "nYriInd = ", nYriInd, "nYriSnp = ", nYriSnp, "\n")
	
	nInd = min(nCeuInd, nYriInd)
	nSnp = min(nCeuInd, nYriSnp)

	nInd = 100
	nSnp = 50
	
	cat("select ", nInd, " individuals, ", nSnp, " snps\n")
		#select 
		ceuSample = ceuGenotype[1:nInd, 1:(2*nSnp)]
		yriSample = yriGenotype[1:nInd, 1:(2*nSnp)]
		
		if(!loadSave)
		{
			cat("calculate ceuSample\n")
			ceuCounts = calculateCounts(ceuSample)
			save(ceuCounts, file = "ceuCounts")		
			cat("calculate yriSample\n")
			yriCounts = calculateCounts(yriSample)
			save(yriCounts, file = "yriCounts")
		}
		else
		{
			load("ceuCounts")
			load("yriCounts")
		}
		
		upperIndex = upper.tri(ceuCounts$c00)
		
		ceuc = ceuCounts
		yric = yriCounts
		
		
		#DEBUG
		cat("ceu c00 dim = ", dim(ceuCounts$c00), "\n")
		
		signd00 = ceuCounts$c00 - yriCounts$c00
		hist(signd00)
		dev.print(device=pdf, "signed_d00.pdf")
		
		#abs difference
		d00 = (abs(ceuCounts$c00 - yriCounts$c00))
		d01 = (abs(ceuCounts$c01 - yriCounts$c01))
		d02 = (abs(ceuCounts$c02 - yriCounts$c02))
		d10 = (abs(ceuCounts$c10 - yriCounts$c10))
		d11 = (abs(ceuCounts$c11 - yriCounts$c11))
		d12 = (abs(ceuCounts$c12 - yriCounts$c12))
		d20 = (abs(ceuCounts$c20 - yriCounts$c20))
		d21 = (abs(ceuCounts$c21 - yriCounts$c21))
		d22 = (abs(ceuCounts$c22 - yriCounts$c22))
		
		#DEBUG
		write.table(ceuCounts$c00, file = "ceuc00")
		write.table(yriCounts$c00, file = "yric00")
		
		hist(ceuCounts$c00)
		dev.print(device=pdf, "ceuc00.pdf")
		hist(d00)
		dev.print(device=pdf, "absd00.pdf")
		hist(yriCounts$c00)
		dev.print(device=pdf, "yric00.pdf")
		
		dev.off()
		
		#DEBUG
		cat("dim(d00)", dim(d00), "\n")
		

		#difference percentage
		p00 = d00/(ceuCounts$c00)
		p01 = d01/(ceuCounts$c01)
		p02 = d02/(ceuCounts$c02)
		p10 = d10/(ceuCounts$c10)
		p11 = d11/(ceuCounts$c11)
		p12 = d12/(ceuCounts$c12)
		p20 = d20/(ceuCounts$c20)
		p21 = d21/(ceuCounts$c21)
		p22 = d22/(ceuCounts$c22)
		
		#DEBUG
		cat("dim(p00)", dim(p00), "\n")
		
		cat("d00 >= c00 ", length(p00[p00>1])/length(p00), "\n")

		
		diff = list("d00" = d00, "d01" = d01, "d02" = d02,
				"d10" = d10, "d10" = d11, "d02" = d12,
				"d20" = d20, "d21" = d21, "d22" = d22)
		
		#print diff
		meandiff = matrix(c(mean(d00, na.rm = T), mean(d01, na.rm = T), mean(d02, na.rm = T), 
						mean(d10, na.rm = T), mean(d11, na.rm = T), mean(d12, na.rm = T),
						mean(d20, na.rm = T), mean(d21, na.rm = T), mean(d22, na.rm = T)), 3, 3, byrow = TRUE)
		
		maxdiff = matrix(c(max(d00, na.rm = T), max(d01, na.rm = T), max(d02, na.rm = T), 
				max(d10, na.rm = T), max(d11, na.rm = T), max(d12, na.rm = T),
				max(d20, na.rm = T), max(d21, na.rm = T), max(d22, na.rm = T)), 3, 3, byrow = TRUE)

		mindiff = matrix(c(min(d00, na.rm = T), min(d01, na.rm = T), min(d02, na.rm = T), 
				min(d10, na.rm = T), min(d11, na.rm = T), min(d12, na.rm = T),
				min(d20, na.rm = T), min(d21, na.rm = T), min(d22, na.rm = T)), 3, 3, byrow = TRUE)

		meanCounts = matrix(c(mean(ceuCounts$c00, na.rm = T), mean(ceuCounts$c01, na.rm = T), mean(ceuCounts$c02, , na.rm = T),
						mean(ceuCounts$c10, na.rm = T), mean(ceuCounts$c11, na.rm = T), mean(ceuCounts$c12, na.rm = T),
						mean(ceuCounts$c20, na.rm = T), mean(ceuCounts$c21, na.rm = T), mean(ceuCounts$c22, na.rm = T)), 3, 3, byrow = T)
		
		
		maxCounts =  matrix(c(max(ceuCounts$c00, na.rm = T), max(ceuCounts$c01, na.rm = T), max(ceuCounts$c02, na.rm = T),
						max(ceuCounts$c10, na.rm = T), max(ceuCounts$c11, na.rm = T), max(ceuCounts$c12, na.rm = T),
						max(ceuCounts$c20, na.rm = T), max(ceuCounts$c21, na.rm = T), max(ceuCounts$c22, na.rm = T)), 3, 3, byrow = T)
		
		pmean = matrix(c(mean(p00, na.rm = T), mean(p01, na.rm = T), mean(p02, na.rm = T), 
						mean(p10, na.rm = T), mean(p11, na.rm = T), mean(p12, na.rm = T),
						mean(p20, na.rm = T), mean(p21, na.rm = T), mean(p22, na.rm = T)), 3, 3)		
		
		
		#pmean = meandiff/meanCounts
		#pmax = maxdiff/maxCounts
		cat("-------------------------------------------------\n")
		cat("counts average\n")
		print(meanCounts)
		cat("counts max\n")
		print(maxCounts)
		
		cat("-------------------------------------------------\n")
		cat("diff\n")
		
		cat("avarage diff\n")
		print(meandiff)
		cat("max diff\n")
		print(maxdiff)

		print(pmean)
		#print(pmax)
}

unitTest <- function()
{
	x = readGenotypeFromFastaFile("/home/xzhou/research_linux/gnome/workspace/data/testdata/x.12encode")
	xc = calculateCounts(x)
	print(xc)
}

testCountsDiff()
#unitTest()

