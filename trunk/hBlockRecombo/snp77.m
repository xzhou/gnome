    cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP';
    clear;
    rawFastaData = fastaread('hapmap_chr7_80SNP_CEU_haplotype.fasta');
    [caseSeq4 refSeq4] = randomSelect(rawFastaData);    
    blocks = [1 15; 16 55; 60 77];
    
    seq = innerBlockDriver(caseSeq4, refSeq4, blocks);
    
    