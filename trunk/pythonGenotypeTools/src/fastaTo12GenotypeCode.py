'''
Created on Jul 14, 2009

@author: xzhou
'''
import re
from encodeFasta import encodeFasta


#preprocess of iregullar file for standard fasta file
def fastaToGenotype(fastaFileName, temp):
    fastaFile = open(fastaFileName)
    lines = fastaFile.readlines()
    
    lines = [line.strip() for line in lines]
    
    nLines = len(lines)
    
    #build RE
    p = re.compile("NA[0-9]*_[AB]")
    
    seq = list()
    
    i = 0
    
    while i < nLines:
        line = ""
        line = lines[i]
        
        if line.startswith(">") and p.match(line[1:len(line)]) != None:
            #line 1 is A
            i = i + 1   #next line
            A = lines[i]
            i = i + 1
            line = lines[i]
            if line.startswith(">") and p.match(line[1:len(line)]) != None:
                i = i + 1
                B = lines[i]
                seq.append(A)
                seq.append(B)
                i = i + 1
            else:
                print "unmatched sequence"
                i = i + 1
        else:
            #continue reading
            i = i + 1
    
    tempFile = open(temp, 'w')
    for i in range(0, len(seq)):
        tempFile.write(seq[i]+"\n")
    tempFile.close()
    #convert genotype sequence
    print i 
    #convert to a genotype file
    
if __name__ == '__main__':
    ceuFileName = "/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/chr10_FGFR2_200kb_phased_CEU.fasta"
    yriFileName = "/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/chr10_FGFR2_200kb_phased_yri.fasta"
    
    ceu = fastaToGenotype(ceuFileName, "/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/ceu.mid")
    encodeFasta("/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/ceu.mid", 
                "/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/ceu.12encode")
    
    yri = fastaToGenotype(yriFileName, "/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/yri.mid")
    encodeFasta("/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/yri.mid", 
                "/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA/yri.12encode")
    
    