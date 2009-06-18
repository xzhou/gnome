'''
Created on Jun 8, 2009

@author: xzhou
'''

'''
convert file to ped data format

combine the random simulated genotype sequence.
'''

def fasta2ped(fastaFileName, pedFileName):
    fastaFile = open(fastaFileName, 'r')
    pedFile = open(pedFileName, 'w')
    
    lines = fastaFile.readlines()
    
    rawSequences = [line.strip() for line in lines if not line.startswith('>')]
    fastaFile.close()
    sequenceLen = len(rawSequences[0])
    nRawSequences = len(rawSequences)
     
    pedSequence = []
    for i in range(0, nRawSequences, 2):
        seq1 = rawSequences[i]
        seq2 = rawSequences[i+1]
        combinedSeq = []
        for j in range(0, len(seq1)):
            combinedSeq.append(seq1[j])
            combinedSeq.append(seq2[j])
        familyInfo = [i, i, 0, 0, 1, 1]     #generating information
        combinedSeq = familyInfo + combinedSeq
        pedSequence.append(combinedSeq)
        
    #output
    for i in range(0, len(pedSequence)):
        for j in range(0, len(pedSequence[i])):
            pedFile.write(str(pedSequence[i][j]))
            pedFile.write(" ")
        pedFile.write("\n")
    pedSequence
    pedFile.close()
    
if __name__ == '__main__':
    fastaFileName = "../../data/sim_4000seq/80SNP_CEU_sim_4000seq.fasta"
    pedFileName = "../../data/sim_4000seq/80SNP_CEU_sim_4000seq.ped"
    fasta2ped(fastaFileName, pedFileName)
    pass