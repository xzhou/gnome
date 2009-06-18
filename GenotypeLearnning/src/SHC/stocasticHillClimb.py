'''
Created on Jun 13, 2009

@author: xzhou
'''

from genotype import *
from math import fabs
import cProfile
from numpy import *
from realRFromFastaFile import readGenotypeFromFasta
import psyco

psyco.full()

def readRFromFile(fileName):
    file = open(fileName, 'r')
    lines = file.readlines()
    
    maxj = 0
    for line in lines:
        line = line.split()
        i = int(line[0])
        j = int(line[1])
        #aRValue = float(line[2])
        if maxj < j:
            maxj = j
    
    realRValues = array([[0.0]*maxj] * maxj)
    #realRValues = zeros((maxj, maxj))
    
    for line in lines:
        line = line.split()
        i = int(line[0])
        j = int(line[1])
        aRValue = float(line[2])
        realRValues[i-1,j-1] = aRValue
        realRValues[j-1,i-1] = aRValue
    
    return realRValues

class SHC(object):
    '''
    SHC is a stochastic hill climbing algorithm
    '''
    
    def __init__(self, max_it, targetGenotype):
        '''
        Constructor
        '''
        self.max_it = 1000  #the max number of iteration
        self.targetGenotype = targetGenotype
        
        self.nSnps = targetGenotype.getNSnps()
        self.nIndividuals = targetGenotype.getNIndividuals()
        self.targetRValue = targetGenotype.calcR(None)
        
        
    def dist(self, targetGenotype, learnnedGenotype):
        rawTarget = targetGenotype.genotype
        rawLearnned = learnnedGenotype.genotype
        
        m, n = shape(rawTarget)
        
        correct = 0.0
        
        for i in range(0, m):
            for j in range(0, n):
                if rawTarget[i,j] == rawLearnned[i,j]:
                    correct = correct + 1
        
        return correct/(m*n)
    
    
    def evaluate(self, sampleRValue):
        '''
        evaluate sample x
        @param x: is a sample
        @return: the fitness of sample x
        '''
        diff = self.targetRValue*self.targetRValue - sampleRValue*sampleRValue
        
        (m,n) = shape(diff)
        
        totalDiff = 0.0
        for i in range(0, m):
            for j in range(0, n):
                aDiff = fabs(diff[i,j])
                totalDiff = totalDiff + aDiff
        return totalDiff

    def singRecoverRate(self, rValues):
        totalSigns = 0
        correctSigns = 0
        (m,n) = shape(rValues)
        if m != n:
            print "invalid r value"
            exit(-1)
        for i in range(0, m):
            for j in range(i+1, n):
                totalSigns = totalSigns + 1
                if(rValues[i,j]*self.targetRValue[i,j] >= 0):
                    correctSigns = correctSigns + 1
        #print correctSigns, "/", totalSigns
        return correctSigns*1.0/totalSigns
    
    def hc(self, max_it, g):
        '''
        hill climbing 5144bd1fd5124e26a7cf4bbe34457cc5246d5178c392eb87b0f6f87b67b66612
        '''
        T = 1
        x = makeGenotype(self.nSnps, self.nIndividuals)
        #quality = self.evaluate(x.calcR(None))
        newDiff = 1
        t = 1
        while t < self.max_it and newDiff != 0:
            print t,
            #generate new sample close to x
            newx = x.singlePointMutation()
            new_r = newx.calcR(None)
            old_r = x.calcR(None)
            newDiff = self.evaluate(new_r)
            oldDiff = self.evaluate(old_r)
            dist = self.dist(self.targetGenotype, newx)
            signRecovered = self.singRecoverRate(new_r)  
            if newDiff < oldDiff:
                x = newx
                print "newDiff = ", newDiff, "oldDiff = ", oldDiff, "signRecovered = ", signRecovered*100, "%", "h_distance=", dist,"%"
                #print newQuality, previousQuality, signRecovered
            t = t+1
    
    def shc(self, max_it, g):
        '''
        statistic hill climbing
        '''
        T = 1
        x = makeGenotype(self.nSnps, self.nIndividuals)
        #quality = self.evaluate(x.calcR(None))
        newQuality = 1
        t = 1
        while t < self.max_it and newQuality != 0:
            print t,
            #generate new sample close to x
            newx = x.singlePointMutation()
            new_r = newx.calcR(None)
            old_r = x.calcR(None)
            newQuality = self.evaluate(new_r)
            oldQuality = self.evaluate(old_r)
            signRecovered = self.singRecoverRate(new_r)
            dist = self.dist(self.targetGenotype, newx)
            p = 1/(1+math.exp((newQuality - oldQuality)/T))
            print "newScore = ", newQuality, "oldQuality = ", oldQuality, "signRecovered = ", signRecovered*100, "%", "h_distance=", dist,"%"
            if random.random() < p :
                x = newx
                #print newQuality, previousQuality, signRecovered
            t = t+1

if __name__ == "__main__":
    #realRValues = readRFromFile("../../data/sim_4000seq/80SNP_CEU_sim_4000seq.rValue");
    psyco.full()
    nSnps = 5
    nIndividuals = 5
    max_it = 1000
    #realRValues = calcualteRFromFasta("../../data/sim_4000seq/80SNP_CEU_sim_4000seq.rValue")
    realGenotype = readGenotypeFromFasta("../../data/sim_4000seq/80SNP_CEU_sim_4000seq.fasta", nSnps, nIndividuals)
    #print realGenotype
    aSHC = SHC(max_it, realGenotype)
    aSHC.hc(aSHC.max_it, 0)
    #cProfile.run(aSHC.shc(aSHC.max_it, 0))
    