%in this file, we will check why I can not calculate r value correctly in
%matlab
cd '/home/xzhou/research_linux/gnome/workspace/data/HAPMAP';

%=========configuration begin====================
fastaFile = 'hapmap_chr7_80SNP_CEU_haplotype.fasta';
%=========configuration end======================


%calculat real R
all = readSeq4(fastaFile);
realR = calcR(all);
alleleMapping = getMajorAllele(all);
[singleAlleleFrequency c] = GnomeCalculator.getSingleAlleleFreq(all, alleleMapping);
singleAlleleFrequency = singleAlleleFrequency';

%estimate R from genotype
genotypeAll = genotypeHelpFuncs.readGenotypeFromFasta(fastaFile);
singleStat = getCounts(genotypeAll);
[totalR pA counts] = estimateR(genotypeAll);
%r from R 
display 'end'