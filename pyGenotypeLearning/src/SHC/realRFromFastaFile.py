'''
Created on Jun 16, 2009

@author: xzhou
'''

from numpy import *
import rpy2.robjects as robjects
import psyco
from genotype import Genotype

psyco.full()
r = robjects.r
r.library("rcalc")


def convert2StandardPEDFormat(rawGenotype):
    nSnps = len(rawGenotype)
    nIndividual = len(rawGenotype[0])/2
    standardGenotype = zeros((nIndividual, 2*nSnps))
    
    for i in range(0, 2*nIndividual, 2):
        for j in range(0, nSnps):
            a = rawGenotype[j][i]
            b = rawGenotype[j][i+1]
            standardGenotype[i/2][2*j] = int(a)
            standardGenotype[i/2][2*j+1] = int(b)
    return standardGenotype

def readGenotypeFromFasta(FileName, nSnps, nIndividuals):
    '''
    @param FileName:     The file name of the fasta File
    @param nSnps:        The limit of snps 
    @param nIndividuals: The limit of individuals
    '''
    file = open(FileName, 'r')
    
    lines = file.readlines()
    
    lines = [line.strip() for line in lines if not line.startswith('>')]
    
    file.close()
    
    realNSnps = len(lines[0])
    realNSeq = len(lines)
    
    if nSnps == -1:
        nSnps = realNSnps
    if nIndividuals == -1:
        nIndividuals = realNSeq/2
    
    if nSnps > realNSnps or 2*nIndividuals > realNSeq:
        print "not enough snps or individuals"
        raise Exception
    
    print nSnps, nIndividuals
    
    lines = [line[0:nSnps] for line in lines]
    
    #rValues = zeros((nSnps, nSnps))
    rValues = array([[0.0]*nSnps]*nSnps)
    rawGenotype = []
    
    
    #majorize for easy compare
    for j in range(0, nSnps):
        aSnps = []
        for i in range(0, nIndividuals):
            seq1 = lines[2*i]
            seq2 = lines[2*i+1]
            a1 = seq1[j]
            a2 = seq2[j]
            #combined = str(a1) + "/" + str(a2)
            aSnps.append(a1)
            aSnps.append(a2)
        rawGenotype.append(aSnps)
    
    #majorize
    genotype = []
    for snp in rawGenotype:
        As = [A for A in snp if A == 'A']
        Cs = [C for C in snp if C == 'C']
        Ts = [T for T in snp if T == 'T']
        Gs = [G for G in snp if G == 'G']   
        nA = len(As)
        nC = len(Cs)
        nT = len(Ts)
        nG = len(Gs)
        counts = [nA, nC, nT, nG]
        counts.sort()
        majorCount = counts[3]
        minorCount = counts[2]
        major = []
        minor = []
        map = {'A':nA, 'C':nC, 'G':nG, 'T':nT}
        for key in map.keys():
            if map[key] == majorCount:
                major.append(key)
            if map[key] == minorCount:
                minor.append(key)
        if len(major) == 1 and len(minor) == 1:
            majorAllele = major[0]
            minorAllele = minor[0]
        elif len(major) == 2:
            majorAllele = major[0]
            minorAllele = major[1]
        elif len(major) == 1 and len(minor) == 3:
            majorAllele = major[0]
            minorAllele = major[0]
        else:
            print "error"
            exit(-1)
        
        for j in range(0, len(snp)):
            if snp[j] == majorAllele:
                snp[j] = "1"
            elif snp[j] == minorAllele:
                snp[j] = "2"

        aFormatedSnp = []
        for j in range(0, len(snp), 2):
            A1 = snp[j]
            A2 = snp[j+1]
            combined = A1 + "/" + A2
            aFormatedSnp.append(combined)
        genotype.append(aFormatedSnp)     
    #format
    print "format complete"
    
    pedGenotype = convert2StandardPEDFormat(rawGenotype)
    
    realGenotype = Genotype(pedGenotype)
    
    return realGenotype
            
if __name__ == '__main__':
    fastaFileName = "../../data/sim_4000seq/80SNP_CEU_sim_4000seq.fasta"
    pedGenotype= readGenotypeFromFasta(fastaFileName, -1, -1)
    pedGenotype.calcR(None)