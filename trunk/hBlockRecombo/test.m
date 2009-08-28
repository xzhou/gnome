%simpleLearning
cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77'

refFileName = 'SIM_100x77_ctl.fasta';
sampleFileName = 'SIM_100x77_smp.fasta';

debug = 1;

targetSeq = readSeq4(sampleFileName);
refSeq = readSeq4(refFileName);

alleleMapping = getMajorAllele(refSeq);

targetR = calcR(targetSeq, alleleMapping);
refR = calcR(refSeq, alleleMapping);

%use the major allele mapping

% b1 = [34 41];   %block one
% b2 = [42 45];   %block two haplotypeFreq = 

%     [hbtarget, ftarget] = getHaplotype(targetSeq, b1);
%     [hbref, fref] = getHaplotype(refSeq, b1);
%     
%     x = [hbtarget ftarget];
%     y = [hbref fref];
% 
%     x = sortrows(x);
%     y = sortrows(y);

hbRecomb(targetR, refSeq, alleleMapping);
