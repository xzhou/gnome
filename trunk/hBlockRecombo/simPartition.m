% partition for simulation data

cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77'

sampleFileName = 'SIM_100x77_smp.fasta';

sampleSeq = fastaread(sampleFileName);

block = hbPartition(sampleSeq);

x = 1;