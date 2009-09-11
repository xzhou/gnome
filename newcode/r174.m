function test = r174()

cd 'F:\IUBResearch\Projects\Bioinfor\data\sim_4000seq'

data = '174SNP_CEU_sim_4000seq.fasta';

seqS = fastaread(data);

nS = length(seqS);
Len = length(seqS(1).Sequence);

int4S = zeros(nS, Len);

for i = 1:length(seqS)
    int4S(i,:) = nt2int(seqS(i).Sequence) - 1;
end

allele1 = int4S(end,:); %why use the last line as allele 0/1 definition
int2S = (int4S == repmat(allele1,nS,1)) + 0;


all_r_S = corrcoef(int2S);
all_r_S(isnan(all_r_S)) = 0;
plotWithSignSimple(all_r_S);

end


