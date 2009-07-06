# TODO: Add comment
# 
# Author: xzhou
###############################################################################


source("SHC.R")

targetGenotype = as.data.frame(matrix(c(1,2,2,1,1,1,2,2,1,1,1,2,2,1,1,1,2,1,1,2), nrow = 1))
sampleGenotype = as.data.frame(matrix(c(1,2,1,2,1,1,1,1,1,1,1,2,2,1,1,1,1,2,2,1), nrow = 1))
#debug(majorize)
#majorize(targetGenotype)

rValues <- calculateRealR(targetGenotype)
sampleR <- calculateRealR(sampleGenotype)
result <- evaluate(rValues,rValues)

print(result)
