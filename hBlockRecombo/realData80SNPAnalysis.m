%%read data LD analysis

cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP';
real80SNPSeq = readSeq4('hapmap_chr7_80SNP_CEU_haplotype.fasta');

majorAllele = getMajorAllele(real80SNPSeq);

rValues = calcR(real80SNPSeq, majorAllele);
rSquare = rValues.*rValues;

h = plotWithSignSimple(rSquare);
saveas(h, 'real80SNPHAPMAPLD.pdf');
