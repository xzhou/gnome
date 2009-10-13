cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP';
fastaData = readSeq4('chr10_FGFR2_200kb_phased_CEU.fasta');

r = calcR(fastaData);

%p = plotWithSignSimple(r);

%saveas(p, 'ceu174.pdf');

blocks = [1 24; 25 45; 46 111; 112 174];