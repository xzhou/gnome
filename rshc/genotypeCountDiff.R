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


#'
#' calculate the halplotype counts of a genotype, note that the genotype is coded such that the haplotype information is preserved
#' @param genotype The 12 encoded genotype that preserve the 1/0 position

calculateHaplotypeCounts <- function(genotype)
{
	genotype <- majorize(genotype)
	m <- nrow(genotype)
	n <- ncol(genotype)
	
	if(n%%2 != 0)
	{
		warning("genotype should be even")
		stop()
	}
	
	nIndividual <- m
	nSnp <- round(n/2)
	
	c00 <- matrix(0.0, nSnp, nSnp)
	diag(c00) <- NA
	c00 <- c01 <- c10 <- c11 <- c00
	
	for(i in 1:(nSnp - 1))
	{
		for(j in (i+1):nSnp)
		{
			hi <- 2*i
			li <- hi - 1
			
			hj <- 2*j
			lj <- hj - 1
			
			A1 <- genotype[,li] #seq 1
			A2 <- genotype[,hi] #seq 2
			
			B1 <- genotype[,lj] #seq 1
			B2 <- genotype[,hj] #seq 2
			
			for(k in 1:nIndividual)
			{
				a1 <- A1[k]
				a2 <- A2[k]
				b1 <- B1[k]
				b2 <- B2[k]
				
				if(a1 == 1 && b1 == 1)
					c00[i,j] = c00[i,j] + 1
				if(a1 == 1 && b1 == 2)
					c01[i,j] = c01[i,j] + 1
				if(a1 == 2 && b1 == 1)
					c10[i,j] = c10[i,j] + 1
				if(a1 == 2 && b1 == 2)
					c11[i,j] = c11[i,j] + 1
				
				
				if(a2 == 1 && b2 == 1)
					c00[i,j] = c00[i,j] + 1
				if(a2 == 1 && b2 == 2)
					c01[i,j] = c01[i,j] + 1
				if(a2 == 2 && b2 == 1)
					c10[i,j] = c10[i,j] + 1
				if(a2 == 2 && b2 == 2)
					c11[i,j] = c11[i,j] + 1
			}	
		}
	}
	
	for(i in 1:nSnp-1)
	{
		for(j in (i+1):nSnp)
		{
			c00[j,i] = c00[i,j]
			c10[j,i] = c10[i,j]
			c01[j,i] = c01[i,j]
			c11[j,i] = c11[i,j]
		}
	}
	
	ret = list("c00" = c00, "c01" = c01, "c10" = c10, "c11" = c11)
	
	#haplotype counts
	class(ret) <- "HC"

	ret
}

testCountsDiff <- function(loadSave = FALSE)
{
	if(!loadSave)
	{
		cat("reading genotype...")
		ceuGenotype = 
				readGenotypeFromFastaFile("/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/ceu.12encode")
		yriGenotype = 
				readGenotypeFromFastaFile("/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/yri.12encode")
		cat("complete\n")
		
		nCeuInd = nrow(ceuGenotype)
		nCeuSnp = ncol(ceuGenotype)/2
		
		nYriInd = nrow(yriGenotype)
		nYriSnp = ncol(yriGenotype)/2
		
		cat("nCeuInd = ", nCeuInd, "nCeuSnp = ", nCeuSnp, "nYriInd = ", nYriInd, "nYriSnp = ", nYriSnp, "\n")
		
		nInd = min(nCeuInd, nYriInd)
		nSnp = min(nCeuInd, nYriSnp)
		
		ceuSample = ceuGenotype
		yriSample = yriGenotype
		
		ceur <- calculateRealR(ceuGenotype)
		yrir <- calculateRealR(yriGenotype)
		
		
		ceuCounts = calculateCounts(ceuSample)
		yriCounts = calculateCounts(yriSample)
		
		ceuHCounts <- calculateHaplotypeCounts(ceuSample)
		yriHCounts <- calculateHaplotypeCounts(yriSample)
		
		
		save(ceuGenotype, file = "ceuGenotype")
		save(yriGenotype, file = "yriGenotype")
		save(ceur, file = "ceur")
		save(yrir, file = "yrir")
		save(ceuCounts, file = "ceuCounts")	
		save(yriCounts, file = "yriCounts")
		save(ceuHCounts, file = "ceuHCounts")
		save(yriHCounts, file = "yriHCounts")
	}
	else
	{
		load(file = "ceuGenotype")
		load(file = "yriGenotype")
		load(file = "ceur")
		load(file = "yrir")
		load(file = "ceuCounts")
		load(file = "yriCounts")
		
		load(file = "ceuHCounts")
		load(file = "yriHCounts")
	}

	
	write.table(ceur, file = "ceur.txt")
	write.table(yrir, file = "yrir.txt")
	
	yris = yrir^2
	ceus = ceur^2
	
	yrisL = yris[yris>0.1]
	ceusL = ceus[ceus>0.1]
	
	cat(">0.01 ceu length = ", length(yrisL)/length(ceus), "yri length = ", length(ceusL)/length(yris))
	
	rsDiff = yris - ceus
	
	cat("\n")
	print(dim(rsDiff))
	cat("\n")
	
	hist(ceur)
	dev.print(device=pdf, "ceur.pdf")
	
	hist(yrir)
	dev.print(device = pdf, "yrir.pdf")
	
	
	hist(rsDiff)
	dev.print(device = pdf, "rsquareDiff.pdf")
	
	plot(ceuHCounts$c00, ceus)
	dev.print(device = pdf, "Hc00vsRS.pdf")
	
	plot(ceuCounts$c00, ceur)
	dev.print(device = pdf, "ceu_c00vsR.pdf")
	
	plot(ceuHCounts$c11, ceus)
	dev.print(device = pdf, "Hc11vsRS.pdf")
	
	plot(abs(ceuHCounts$c00*ceuHCounts$c11), ceur)
	dev.print(device = pdf, "abshc0011vsR.pdf")
	
	rsdiff = ceus - yris
	cdiff = ceuCounts$c00 - yriCounts$c00
	plot(cdiff, rsdiff)
	dev.print(device=pdf, file = "rdiff_cdiff.pdf")
	
	rdiff = ceur - yrir
	c00_x_c11diff = ceuCounts$c00*ceuCounts$c11 - yriCounts$c00*yriCounts$c11
	plot(c00_x_c11diff, rsdiff)
	dev.print(device=pdf, file = "rsdiff_c0011diff.pdf")
	
	L = (ceuHCounts$c00+ceuHCounts$c01)*(ceuHCounts$c11+ceuHCounts$c10)*
			(ceuHCounts$c00+ceuHCounts$c10)*(ceuHCounts$c11+ceuHCounts$c01)
	
	plot((ceuHCounts$c00*ceuHCounts$c11-ceuHCounts$c01*ceuHCounts$c10)/sqrt(L), ceus)
	dev.print(device = pdf, "rVSrs.pdf")
	
	c001122 = ceuCounts$c00 + ceuCounts$c01/2 + ceuCounts$c22
	plot(c001122, ceur)
	dev.print(device = pdf, "c001122vsR.pdf")

		upperIndex = upper.tri(ceuCounts$c00)
		
		ceuc = ceuCounts
		yric = yriCounts
		
		
		#DEBUG
		cat("ceu c00 dim = ", dim(ceuCounts$c00), "\n")
		
		signd00 = ceuCounts$c00 - yriCounts$c00
		hist(signd00)
		dev.print(device=pdf, "signed_d00.pdf")
		
		#abs difference
		d00 = ceuCounts$c00 - yriCounts$c00
		d01 = ceuCounts$c01 - yriCounts$c01
		d02 = ceuCounts$c02 - yriCounts$c02
		d10 = ceuCounts$c10 - yriCounts$c10
		d11 = ceuCounts$c11 - yriCounts$c11
		d12 = ceuCounts$c12 - yriCounts$c12
		d20 = ceuCounts$c20 - yriCounts$c20
		d21 = ceuCounts$c21 - yriCounts$c21
		d22 = ceuCounts$c22 - yriCounts$c22
		
		#DEBUG
		write.table(ceuCounts$c00, file = "ceuc00")
		write.table(yriCounts$c00, file = "yric00")
		
		hist(ceuCounts$c00)
		dev.print(device=pdf, "ceuc00.pdf")
		hist(ceuCounts$c01)
		dev.print(device=pdf, "ceuc01.pdf")
		hist(ceuCounts$c02)
		dev.print(device=pdf, "ceuc02.pdf")	
		hist(ceuCounts$c10)
		dev.print(device=pdf, "ceuc10.pdf")		
		hist(ceuCounts$c11)
		dev.print(device=pdf, "ceuc11.pdf")		
		hist(ceuCounts$c12)
		dev.print(device=pdf, "ceuc12.pdf")		
		hist(ceuCounts$c20)
		dev.print(device=pdf, "ceuc20.pdf")	
		hist(ceuCounts$c21)
		dev.print(device=pdf, "ceuc21.pdf")
		hist(ceuCounts$c22)
		dev.print(device=pdf, "ceuc22.pdf")
		
		hist(yriCounts$c00)
		dev.print(device=pdf, "yric00.pdf")
		
		hist(d00)
		dev.print(device=pdf, "d00.pdf")
		hist(d01)
		dev.print(device=pdf, "d01.pdf")
		hist(d02)
		dev.print(device=pdf, "d02.pdf")
		hist(d10)
		dev.print(device=pdf, "d10.pdf")
		hist(d11)
		dev.print(device=pdf, "d11.pdf")
		hist(d12)
		dev.print(device=pdf, "d12.pdf")
		hist(d20)
		dev.print(device=pdf, "d20.pdf")
		hist(d21)
		dev.print(device=pdf, "d21.pdf")
		hist(d22)
		dev.print(device=pdf, "d22.pdf")

		hist(yriCounts$c00)
		dev.print(device=pdf, "yric00.pdf")
		
		dev.off()
		
		#DEBUG
		cat("dim(d00)", dim(d00), "\n")
		index <- abs(d00) > 30
		large_c00_rs_diff = rsDiff[index]
		hist(large_c00_rs_diff)
		dev.print(device = pdf, "large_c00_rs_diff.pdf")
		
		large_c00_ceu_rs = ceus[index]
		hist(large_c00_ceu_rs)
		dev.print(device=pdf, "large_c00_ceu_rs.pdf")
		
		large_c00_yri_rs = yris[index]
		hist(large_c00_yri_rs)
		dev.print(device = pdf, "large_c00_yri_rs.pdf")
		
		large_c00_ceu_c00 = (ceuCounts$c00)[index]
		hist(large_c00_ceu_c00)
		dev.print(device = pdf, "large_c00_ceu_c00.pdf")
		
		large_c00_ceu_c22 = (ceuCounts$c22)[index]
		dev.print(device = pdf, "large_c00_ceu_c22.pdf")
		
		large_c00_yri_c00 = (yriCounts$c00)[index]
		hist(large_c00_yri_c00)
		dev.print(device = pdf, "large_c00_yri_c00.pdf")
		
		plot(large_c00_yri_c00, large_c00_yri_rs)
		dev.print(device = pdf, "large_c00_yri_c00_rs.pdf")
		
		plot(large_c00_ceu_c00, ceus[index])
		dev.print(device = pdf, "large_c00_ceu_c00_rs.pdf")
		
		#difference percentge
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
	
	print(x)
	
	xc = calculateCounts(x)
	hc <- calculateHaplotypeCounts(x)
	r = calculateRealR(x)
	#print(xc)
	#print(r)
	cat("\n---------------\n")
	print(hc)
}

testCountsDiff(T)
#unitTest()

