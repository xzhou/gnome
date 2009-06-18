'''
Created on Jun 13, 2009

@author: xzhou
'''

from numpy import *     #import matrix support
from rpy2 import *      #import R interface and calling optimize function
import rpy2.robjects as robjects
import copy


def fillGenotypeWithRandomValue(genotype):
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

def makeGenotype(nSnps, nIndividual):
    nSnps = nSnps;
    nIndividual = nIndividual
    genotype = array([[0]*(2*nSnps)]*nIndividual)   #we use the ped format
    fillGenotypeWithRandomValue(genotype)
    return Genotype(genotype)




class Genotype(object):
    '''
    GenotypeSample represent a sample of binary encoded genotype sequence
    ''' 
    r = robjects.r
    r.library("rcalc")
    
    def __init__(self, genotype):
        (m,n) = shape(genotype)
        self.nSnps = n/2
        self.nIndividual = m
        self.genotype = genotype

    def getNSnps(self):
        return self.nSnps
    
    def getNIndividuals(self):
        return self.nIndividual
    
    def setGenotype(self, genotype):
        self.genotype = genotype
        (m,n) = shape(genotype)
        self.nSnps = n/2
        self.nIndividual = m
        
    def convert2RGenotypeFormat(self, snp1):
        newsnp1 = []
        for i in range(0, len(snp1[0])):
            locus1 = str(snp1[0][i]) + "/" + str(snp1[1][i])
            newsnp1.append(locus1)
        return newsnp1

    def calcR2Snps(self, snp1, snp2):
        g1 = Genotype.r.genotype(snp1)
        g2 = Genotype.r.genotype(snp2)
        result = Genotype.r.rcalc(g1, g2)
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
        
        
        
        rValues = array([[2.0]*nSnps]*nSnps);
        for i in range(0, len(formatedSnps)-1):
            for j in range(i+1, len(formatedSnps)):
                snp1 = formatedSnps[i]
                snp2 = formatedSnps[j]
                r = self.calcR2Snps(snp1, snp2)
                rValues[i,j] = r
                rValues[j,i] = r
                #print i,j,r
        return rValues

    def singlePointMutation(self):
        
        genotype = self.genotype
        #create a new copy
        newSample = copy.deepcopy(self)
        
        x = random.randint(0, newSample.nIndividual)
        y = random.randint(0, 2*newSample.nSnps)
        
        originalValue = newSample.genotype[x,y]
        if originalValue == 1:
            newSample.genotype[x,y] = 2
        else:
            newSample.genotype[x,y] = 1
        
        return newSample
    
    def segmentMutation(self, genotype):
        if genotype == None:
            genotype = self.genotype

if __name__ == "__main__":
    aSample = Genotype(77, 100)
    rValues = aSample.calcR(aSample.genotype)
    
    