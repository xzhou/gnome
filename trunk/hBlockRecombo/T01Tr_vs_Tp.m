%this script will test the power of Tr and Homer's test over the simulated
%data

cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77';

caseSeq = fastaread('80SNP_CEU_sim_4000seq_control1000.fasta');


T = 0;

powerCurve(caseSeq, T, 200);
