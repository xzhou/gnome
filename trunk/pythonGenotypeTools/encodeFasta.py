'''
Created on Jun 19, 2009

@author: xzhou
'''
import psyco
psyco.full()

def encode0112(fileName, outputFileName, hapFileName = "", nSnps = -1, nIndividuals = -1):
    file = open(fileName, 'r')
    
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
    
    rawGenotype = []

    #majorize for easy compare
    for j in range(0, nSnps):
        aSnps = []
        for i in range(0, nIndividuals):
            seq1 = lines[2*i]
            seq2 = lines[2*i+1]
            a1 = seq1[j]
            a2 = seq2[j]
            aSnps.append(a1)
            aSnps.append(a2)
        rawGenotype.append(aSnps)
    #print rawGenotype[0]
    
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
                snp[j] = 1
            elif snp[j] == minorAllele:
                snp[j] = 2
        genotype.append(snp)
    
    if hapFileName != "":
        hapFile = open(hapFileName, 'w')
        
        for i in range(0, 2*nIndividuals):
            for j in range(0, len(genotype)):
                a = genotype[j][i]
                hapFile.write(str(a-1) + " ")
                
            hapFile.write("\n");

    outputFile = open(outputFileName, 'w')
    for j in range(0, 2*nIndividuals, 2):
        for i in range(0, len(genotype)):
            a1 = genotype[i][j]
            a2 = genotype[i][j+1]
            outputFile.write(str(a1) + " " + str(a2) + " ")
        outputFile.write("\n")
    outputFile.close()


if __name__ == '__main__':
    fastaFileName = "../../data/sim_4000seq/80SNP_CEU_sim_4000seq.fasta"
    outputFileName = "../../data/sim_4000seq/80SNP_CEU_sim_4000seq.12encode"
    hapFileName = "../../data/sim_4000seq/80SNP_CEU_sim_4000seq.01hap"
    encode0112(fastaFileName, outputFileName, hapFileName)
    