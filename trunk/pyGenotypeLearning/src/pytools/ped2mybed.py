'''
Created on Jun 7, 2009

@author: xzhou
'''

def convertped2mybed(pedFileName, mybedFileName):
    pedFile = open(pedFileName);
    
    lines = pedFile.readlines();
    
    individualSequences = []    #the genotype sequcence of a individual
    
    '''read the genotype data'''
    for line in lines:
        tokens = line.split()
        genotype = tokens[6:len(tokens)-1]
        individualSequences.append(genotype)
    
    aGenotype = individualSequence[0]
    print len(aGenotype)
    
    snps = []
    
    #convert to snps
    for i in range(0, len(aGenotype)):
        aSnp = [];
        for genotype in individualSequences:
            aSnp.append(genotype[i])
        snps.append(aSnp)
    
    for i in range(0, len(snps), 2):
        snp1 = snps[i]
        snp2 = snps[i+1]
        
        length = len(snp1) #assume of the same length
        
        As = [A for A in snp1 if A == 'A'] + [A for A in snp2 if A == 'A'];
            
                   
    
        
    
    
    
    

if __name__ == '__main__':
    
    fileName = "../data/77SNPHapmap.ped"
    mybedFileName = "../data/77SNPHapmap.mybed"
    
    convertped2mybed(pedFileName, mybedFileName);