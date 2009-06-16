'''
Created on Jun 13, 2009

@author: xzhou
'''

from numpy import *     #import matrix support
from rpy2 import *      #import R interface and calling optimize function
import rpy2.robjects as robjects
import copy

class GenotypeSample(object):
    '''
    GenotypeSample represent a sample of binary encoded genotype sequence
    ''' 
    r = robjects.r
    r.library("rcalc")
    
    def __init__(self, nSnps, nIndividual):
        self.nSnps = nSnps;
        self.nIndividual = nIndividual
        self.genotype = array([[0]*(2*nSnps)]*nIndividual)   #each snp has two allele
        self.generateRandomGenotypeSample(self.genotype)

    
    def setGenotype(self, genotype):
        self.genotype = genotype
        
    def convert2RGenotypeFormat(self, snp1):
        newsnp1 = []
        for i in range(0, len(snp1[0])):
            locus1 = str(snp1[0][i]) + "/" + str(snp1[1][i])
            newsnp1.append(locus1)
        return newsnp1

    def calcR2Snps(self, snp1, snp2):
        g1 = GenotypeSample.r.genotype(snp1)
        g2 = GenotypeSample.r.genotype(snp2)
        result = GenotypeSample.r.rcalc(g1, g2)
        rValue = result[0][0]
        #print "r=", result[0][0], "pAB=", result[1][0]
        return rValue
    
    #convert to R format like ("A/B", "B/B")
    def formateSnps(self):
        (m,n) = shape(self.genotype)
        formatedSnps = []
        for i in range(0, n, 2):
            snp1 = [self.genotype[:,i], self.genotype[:, i+1]]
            formatedSnp = self.convert2RGenotypeFormat(snp1)
            formatedSnps.append(formatedSnp)
        return formatedSnps

    def calcR(self, genotype):
        '''
            calculate the r value of this sample
        '''
        if genotype == None:
            genotype = self.genotype 
        (m,n) = shape(genotype)
        #test if n is a even number since n is 2*nSnps
        if n%2 != 0:
            return
        nSnps = n/2
        formatedSnps = self.formateSnps()
        rValues = array([[2]*nSnps]*nSnps);
        for i in range(0, len(formatedSnps)-1):
            for j in range(i+1, len(formatedSnps)):
                snp1 = formatedSnps[i]
                snp2 = formatedSnps[j]
                r = self.calcR2Snps(snp1, snp2)
                rValues[i,j] = r
                rValues[j,i] = r
        return rValues
    
    def generateRandomGenotypeSample(self, genotype):
        '''
        @param genotype:    the genotype
        '''
        (n,m) = shape(genotype)
        for i in range(0, n):
            for j in range(0,m):
                x = random.random()
                if x <= 0.5:
                    genotype[i,j] = 1       #major
                else:
                    genotype[i,j] = 2       #minor
        return genotype

    def singlePointMutation(self, genotype):
        if genotype == None:
            genotype = self.genotype
            
        '''
        newSample = GenotypeSample(77, 100)
        newSample.genotype = self.genotype
        newSample.nSnps = self.nSnps
        newSample.nIndividual = self.nIndividual
        '''
        #create a new copy
        newSample = copy.deepcopy(self)
        
        x = random.randint(0, newSample.nIndividual)
        y = random.randint(0, newSample(2*nSnps))
        
        originalValue = newSample.genotype[i,j]
        if originalValue == 1:
            newSample.genotype[i,j] = 2
        else:
            newSample.genotype[i,j] = 1
        
        return newSample
    
    def segmentMutation(self, genotype):
        if genotype == None:
            genotype = self.genotype
        
        
    
if __name__ == "__main__":
    aSample = GenotypeSample(77, 100)
    rValues = aSample.calcR(aSample.genotype)
    
    