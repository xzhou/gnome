

blocks = [1 2 3 4 5 6 7 9  12 18 22 23 24 30 31 32 34 42 46 50 51 53 54 58 61 65 68 69 71 72 73 74; 
          1 2 3 4 5 6 8 11 17 21 22 23 29 30 31 33 41 45 49 50 52 53 57 60 64 67 68 70 71 72 73 77]';

%blocks = [34 41; 42 45];
      
cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77'

refFileName = 'SIM_100x77_ctl.fasta';
sampleFileName = 'SIM_100x77_smp.fasta';

caseSeq4 = readSeq4(sampleFileName);
refSeq4 = readSeq4(refFileName);

%recursiveLearning(caseSeq4, refSeq4, blocks);
%recursiveLearning(caseSeq4, caseSeq4, blocks);

findMultipleSolution(caseSeq4, caseSeq4, blocks);