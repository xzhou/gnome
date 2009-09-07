%this file test the power of SIM data
cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77'

caseFastaFile = 'SIM_100x77_smp.fasta';
refFastaFile = 'SIM_100x77_ctl.fasta';

seqS = fastaread(caseFastaFile);
seqR = fastaread(refFastaFile);

