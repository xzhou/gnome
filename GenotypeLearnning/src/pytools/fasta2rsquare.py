'''
Created on Jun 10, 2009

@author: xzhou

'''

'''
individuals:     lines of individual haplotype sequence
return value:    the snps
'''

import math

def individual2snps(individuals):
    
    nIndividual = len(individuals)

    nSnps = len(individuals[0])
    
    snps = []
    for j in range(0, nSnps):
        aSnp = []
        for i in range(0, nIndividual):
            aSnp.append(individuals[i][j])
        snps.append(aSnp)
    return snps

def encodeSnps(snps):
    for i in range(0, len(snps)):
        aSnp = snps[i]
        As = [A for A in aSnp if A == 'A']
        Cs = [C for C in aSnp if C == 'C']
        Gs = [G for G in aSnp if G == 'G']
        Ts = [T for T in aSnp if T == 'T']
        
        pureList = [As, Cs, Gs, Ts]
        
        nA = len(As)
        nC = len(Cs)
        nT = len(Ts)
        nG = len(Gs)
        
        counts = [nA, nT, nC, nG]
        counts.sort()
        
        majorCount = counts[3]
        minorCount = counts[2]
        
        major = 'X'
        minor = 'x'
        
        #find major and minor allele
        if majorCount != minorCount:
            nameMap = {nA:'A', nC:'C', nG:'G', nT:'T'}
            major = nameMap[majorCount]
            minor = nameMap[minorCount]
        else:
            tmp = []
            if nA == majorCount:
                tmp.append('A')
            if nC == majorCount:
                tmp.append('C')
            if nT == majorCount:
                tmp.append('T')
            if nG == majorCount:
                tmp.append['G']
            
            major = tmp[0]
            minor = tmp[1]
        
        for j in range(0, len(aSnp)):
            if aSnp[j] == major:
                aSnp[j] = 0
            elif aSnp[j] == minor:
                aSnp[j] = 1
            else:
                raise Exception
        
def calculateRSquares(snps, rFileName, rsFileName):
    r = []
    rs = []
    
    rFile = open(rFileName, 'w')
    rsFile = open(rsFileName, 'w')
    
    nSnps = len(snps)
    for i in range(0, nSnps-1):
        for j in range(i+1, nSnps):
            snp1 = snps[i]
            snp2 = snps[j]
            c00 = 0
            c01 = 0
            c10 = 0
            c11 = 0
            for k in range(0, len(snp1)):
                if snp1[k] == 0 and snp2[k] == 0:
                    c00 = c00 + 1
                elif snp1[k] == 0 and snp2[k] == 1:
                    c01 = c01 + 1
                elif snp1[k] == 1 and snp2[k] == 0:
                    c10 = c10 + 1
                elif snp1[k] == 1 and snp2[k] == 1:
                    c11 = c11 + 1
            c0x = c00 + c01
            c1x = c10 + c11
            cx0 = c00 + c10
            cx1 = c01 + c11
            
            D = c00*c11 - c01*c10
            L = c0x*c1x*cx0*cx1
            
            if L == 0:
                pass
            
            r_ij = D*1.0/math.sqrt(L)
            rs_ij = r_ij * r_ij
            
            rList = [i, j, r_ij]
            rsList = [i, j, rs_ij]
            
            rFile.write(str(i+1) + " " + str(j+1) + " " + str(r_ij) + "\n")
            rsFile.write(str(i+1) + " " + str(j+1) + " " + str(rs_ij) + "\n")
            
            r.append(rList)
            rs.append(rsList)
            #end for k
        #end for j
    #end for i
    
    rFile.close()
    rsFile.close()
    
    return [r, rs]
                

def fasta2rsquare(fastaFileName, rFileName, rsFileName):
    fastaFile = open(fastaFileName, 'r')
    lines = fastaFile.readlines()
    lines = [line.strip() for line in lines if not line.startswith(">")]
    fastaFile.close()
    
    snps = individual2snps(lines)
    encodeSnps(snps)
    
    [r, rs] = calculateRSquares(snps, rFileName, rsFileName)
    
    
    
if __name__ == '__main__':
    fastaFileName = "../../data/sim_4000seq/80SNP_CEU_sim_4000seq.fasta"
    rFileName = "../../data/sim_4000seq/80SNP_CEU_sim_4000seq.rValue"
    rsFileName = "../../data/sim_4000seq/80SNP_CEU_sim_4000seq.rsValue"
    fasta2rsquare(fastaFileName, rFileName, rsFileName)
    pass