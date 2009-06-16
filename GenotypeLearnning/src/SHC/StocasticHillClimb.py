'''
Created on Jun 13, 2009

@author: xzhou
'''

import genotypeSample
import math

class SHC(object):
    '''
    SHC is a stochastic hill climbing algorithm
    '''
    def __init__(self, realRValue):
        '''
        Constructor
        '''
        self.max_it = 1000  #the max number of iteration
        self.realRValue =  realRValue
        
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
        diff = self.realRValue*self.realRValue - sampleRValue*sampleRValue
        
        (m,n) = shape(diff)
        
        totalDiff = 0.0
        for i in range(0, m):
            for j in range(0, n):
                aDiff = math.abs(diff[i,j])
                totalDiff = totalDiff = aDiff
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
                totalSings = totalSigns + 1
                if(rValues[i,j]*self.realRValue[i,j] >= 0):
                    correctSigns = correctSigns + 1
        return correctSigns*1.0/totalSigns
        
    
    def shc(self, max_it, g):
        '''
        statistic hill climbing
        '''
        T = 1
        
        x = genotypeSample(77, 100)
        quality = evalueate(x.calcR)
        t = 1
        while t < self.max_it :
            newx = x.singlePointMutation()
            newQuality = evaluate(newx.calcR())
            previousQuality = evaluate(x.calcR())
            
            p = 1/(1+math.exp((newQuality - previousQuality)/T))
            
            if random.random(0,1) < p :
                x = newx
                print newQuality, previousQuality
            t = t+1
        