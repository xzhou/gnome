%% compare sim and read data

cd '/home/xzhou/research_linux/gnome/workspace/data/SIM_HAPMPA_COMP';

SNP80_SIM_CEU = readSeq4('80SNP_CEU_sim_4000seq.fasta');
SNP80_REAL_CEU = readSeq4('hapmap_chr7_80SNP_CEU_haplotype.fasta');

SNP174_SIM_CEU = readSeq4('174SNP_CEU_sim_4000seq.fasta');
SNP174_REAL_CEU = readSeq4('chr10_FGFR2_200kb_phased_CEU.fasta');

SIM_80_CEU_R = calcR(SNP80_SIM_CEU);
REAL_80_CEU_R = calcR(SNP80_REAL_CEU);

SIM_174_CEU_R = calcR(SNP174_SIM_CEU);
REAL_174_CEU_R = calcR(SNP174_REAL_CEU);

h = plotWithSignSimple(SIM_80_CEU_R, REAL_80_CEU_R, 0, 0);
ylabel('REAL 80 CEU_R');
xlabel('SIM 80 CEU_R');
saveas(h, '80SNP SIM REAL.pdf');

h = plotWithSignSimple(SIM_174_CEU_R, REAL_174_CEU_R, 0, 0);
ylabel('REAL 174 CEU_R');
xlabel('SIM 174 CEU_R');
saveas(h, '174SNP SIM REAL.png');



