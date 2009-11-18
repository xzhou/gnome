'''
Created on Jun 8, 2009
rafdfdf
@author: xzhou
'''


def fasta2genotype(fastaFileName, genotypeFileName):
    fastaFile = open(fastaFileName)
    lines = fastaFile.readlines()
    
    individualSequences = []
    
    '''read the genotype data'''
    for line in lines:
        if line.startswith(">"):
            continue
        individualSequences.append(line.strip())
    
    '''get the length of the sequence, assume all the sequence has the same len'''
    aSampleSequence = individualSequences[0]
    seqLen = len(aSampleSequence)
    
    snps = []
    
    '''extract to column mode'''
    for i in range(0, seqLen):
        aSnp = []
        for seq in individualSequences:
            aSnp.append(seq[i])
        snps.append(aSnp)
    
    
    '''determine major and minor'''
    for i in range(0, len(snps)):
        snp1 = snps[i]  
        snpLen = len(snp1)
        
        As = [A for A in snp1 if A == 'A']
        Ts = [T for T in snp1 if T == 'T']
        Cs = [C for C in snp1 if C == 'C']
        Gs = [G for G in snp1 if G == 'G']
        
        nA = len(As) #count of A
        nT = len(Ts)
        nC = len(Cs)
        nG = len(Gs)
        
        counts = [nA, nT, nC, nG]
        counts.sort()
        
        majorCount = counts[3]
        minorCount = counts[2]
        
        '''in case major and minor are the same'''
        if majorCount == minorCount:
            minorCount = majorCount - 1
        
        nameMap = {nA:'A', nC:'C', nG:'G', nT:'T'}
        
        if len(nameMap) != 3:
            print snp1
        
        #print nameMap, majorCount, minorCount
        
        major = nameMap[majorCount]
        minor = nameMap[minorCount]

        for i in range(0, snpLen):
            if snp1[i] == major:
                snp1[i] = 0
            elif snp1[i] == minor:
                snp1[i] = 1
            else:
                print "error"
    
    snpGenotype = []
    
    '''convert genotype'''
    for i in range(0, len(snps[0]), 2):
        aSnpGenotype = []
        for j in range(0, len(snps)):
            allele1 = snps[j][i]
            allele2 = snps[j][i+1]
            #print allele1, allele2     
            x = allele1 + allele2
            aSnpGenotype.append(x)
        
        snpGenotype.append(aSnpGenotype)
    
    '''print output'''
    outFile = open(genotypeFileName, 'w')
    for i in range(0, len(snpGenotype)):
        for j in range(0, len(snpGenotype[i])):
            outFile.write(str(snpGenotype[i][j]))
        outFile.write("\n")
        
    outFile.close()
    fastaFile.close()
if __name__ == '__main__':
    fastaFileName = "../../data/sim_4000seq/80SNP_CEU_sim_4000seq.fasta"
    genotypeFileName = "../../data/sim_4000seq/80SNP_CEU_sim_4000seq.genotype"
    fasta2genotype(fastaFileName, genotypeFileName)
