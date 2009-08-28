cd '/home/xzhou/research_linux/gnome/workspace/data/88_77_CEU_YRI_DATA';

ceuFileName = 'chr10_FGFR2_200kb_phased_CEU.fasta';
yriFileName = 'chr10_FGFR2_200kb_phased_yri.fasta';

seqCEU = fastaread(ceuFileName);
seqYRI = fastaread(yriFileName);

mCEU = length(seqCEU);
nCEU = length(seqCEU(1).Sequence);

mYRI = length(seqYRI);
nYRI = length(seqYRI(1).Sequence);

int4CEU = zeros(mCEU, nCEU);
int4YRI = zeros(mYRI, nYRI);

parfor i = 1:mCEU
    int4CEU(i,:) = nt2int(seqCEU(i).Sequence) - 1;
end

parfor i = 1:mYRI
    int4YRI(i,:) = nt2int(seqYRI(i).Sequence) - 1;
end


%define allele ???
allele1 = int4CEU(end,:);
ceuMapping = getMajorAllele(int4CEU);
yriMapping = getMajorAllele(int4YRI);

int2CEU = (int4CEU == repmat(ceuMapping, mCEU, 1)) + 0;
int2YRI = (int4YRI == repmat(ceuMapping, mYRI, 1)) + 0;

all_r_CEU = corrcoef(int2CEU);
all_r_CEU(isnan(all_r_CEU)) = 0;

all_r_YRI = corrcoef(int2YRI);
all_r_YRI(isnan(all_r_YRI)) = 0;

%plot the graph
[m n] = size(all_r_YRI);
index = ones(m);
index = logical(triu(index));
mixR = all_r_YRI;
mixR(index) = all_r_CEU(index);

%h = plotWithSignSimple(all_r_CEU, all_r_YRI, 0.0);
%saveas(h, 'ceuyri.pdf');

h = plotWithSignSimple(all_r_CEU(1:80, 1:80), all_r_YRI(1:80, 1:80),  0.1);
saveas(h, 'ceuyri_1_80.pdf');

h = plotWithSignSimple(all_r_CEU(75:end, 75:end), all_r_YRI(75:end, 75:end), 0.1);
saveas(h, 'ceuyri_75_174.pdf');



