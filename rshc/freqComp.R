#function get block frequency
getBlockFreq <- function(genotype)
  {
  }



#This function will test if the genotype frequency of two samples are the same
freqComp <- function(caseGenotype, refGenotype)
  {
    
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
