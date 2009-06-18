'''
Created on Jun 13, 2009

@author: xzhou
'''

from genotypeSample import GenotypeSample
from math import fabs
import cProfile
from numpy import *
from realRFromFastaFile import calcualteRFromFasta
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
    
    def __init__(self, targetGenotype, realRValue, nSnps, nIndividuals):
        '''
        Constructor
        '''
        self.max_it = 1000  #the max number of iteration
        self.realRValue =  realRValue
        self.nSnps = nSnps
        self.nIndividuals = nIndividuals
        self.targetGenotype = targetGenotype
        
        (m,n) = shape(realRValue)
        if m != n or m != nSnps:
            print "error dimensions"
            exit(-1)

    def initialize(self):
        '''
        generate initial sample
        '''
        pass
    def evaluate(self, sampleRValue):
        '''
        evaluate sample x
        @param x: is a sample
        @return: the fitness of sample x
        '''
        #calcuate the difference between r square
        #print shape(self.realRValue), shape(sampleRValue)
        
        diff = self.realRValue*self.realRValue - sampleRValue*sampleRValue
        
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
                if(rValues[i,j]*self.realRValue[i,j] >= 0):
                    correctSigns = correctSigns + 1
        #print correctSigns, "/", totalSigns
        return correctSigns*1.0/totalSigns
    
    def hc(self, max_it, g):
        '''
        hill climbing 5144bd1fd5124e26a7cf4bbe34457cc5246d5178c392eb87b0f6f87b67b66612
        '''
        T = 1
        x = GenotypeSample(self.nSnps, self.nIndividuals)
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
            signRecovered = self.singRecoverRate(new_r)  
            if newDiff < oldDiff:
                x = newx
                print "newDiff = ", newDiff, "oldDiff = ", oldDiff, "signRecovered = ", signRecovered*100, "%"
                #print newQuality, previousQuality, signRecovered
            t = t+1
    
    def shc(self, max_it, g):
        '''
        statistic hill climbing
        '''
        T = 1
        x = GenotypeSample(self.nSnps, self.nIndividuals)
        #quality = self.evaluate(x.calcR(None))
        newQuality = 1
        t = 1
        while t < self.max_it and newQuality != 0:
            print t,
            #generate new sample close to x
            newx = x.singlePointMutation()
            new_r = newx.calcR(None)
            old_r = x.calcR(None)
            dist = 
            newQuality = self.evaluate(new_r)
            oldQuality = self.evaluate(old_r)
            signRecovered = self.singRecoverRate(new_r)
            p = 1/(1+math.exp((newQuality - oldQuality)/T))
            print "newScore = ", newQuality, "oldQuality = ", oldQuality, "signRecovered = ", signRecovered*100, "%"
            if random.random() < p :
                x = newx
                #print newQuality, previousQuality, signRecovered
            t = t+1

if __name__ == "__main__":
    #realRValues = readRFromFile("../../data/sim_4000seq/80SNP_CEU_sim_4000seq.rValue");
    psyco.full()
    nSnps = 5
    nIndividuals = 5
    #realRValues = calcualteRFromFasta("../../data/sim_4000seq/80SNP_CEU_sim_4000seq.rValue")
    realRValues, realGenotype = calcualteRFromFasta("../../data/sim_4000seq/80SNP_CEU_sim_4000seq.fasta", nSnps, nIndividuals)
    #print realGenotype
    aSHC = SHC(realGenotype, realRValues, nSnps, nIndividuals)
    aSHC.hc(aSHC.max_it, 0)
    #cProfile.run(aSHC.shc(aSHC.max_it, 0))
    