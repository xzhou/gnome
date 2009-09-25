function test2 = strongWeakAnalysis()

cd 'D:\IUBResearch\Projects\Bioinfor\data\88_77_CEU_YRI_DATA'

all = 'hapmap_chr7_80SNP_CEU_haplotype.fasta';

seqAll = fastaread(all);

seqAll(end) = [];
seqAll(end) = [];

sel = [1:116];
a = 116;
b = randi(a,1,1);
selS = zeros(1, 58);
selR = zeros(1, 58);

for i = 1:58
    c = randi(a,1,1);
    selS(1,i) = sel(1, c);
    sel(c) = [];
    a = a-1;
end

selR = sel;

for i = 1:58
    seqS(i*2-1) = seqAll(selS(i)*2-1);
    seqS(i*2) = seqAll(selS(i)*2);
    
    seqR(i*2-1) = seqAll(selR(i)*2-1);
    seqR(i*2) = seqAll(selR(i)*2);
end

seqS = seqS';
seqR = seqR';

nS = length(seqS);
Len = length(seqS(1).Sequence);
nR = length(seqR);

int4S = zeros(nS, Len);
int4R = int4S;

for i = 1:length(seqS)
    int4S(i,:) = nt2int(seqS(i).Sequence) - 1;
    int4R(i,:) = nt2int(seqR(i).Sequence) - 1;
end

allele1 = int4S(end,:); %why use the last line as allele 0/1 definition
int2S = (int4S == repmat(allele1,nS,1)) + 0;
int2R = (int4R == repmat(allele1,nR,1)) + 0;

all_r_S = corrcoef(int2S);
all_r_S(isnan(all_r_S)) = 0;

all_r_R = corrcoef(int2R);
all_r_R(isnan(all_r_R)) = 0;

StatS.Tr = zeros(nS,1);
StatR.Tr = zeros(nS,1);

for i = 1:nS
    StatS.Tr(i) = getTr(int2S(i,:), all_r_S, all_r_R);
    StatR.Tr(i) = getTr(int2R(i,:), all_r_S, all_r_R);
end
StatS.Tr = StatS.Tr/sqrt(Len*(Len-1)/2);  %??
StatR.Tr = StatR.Tr/sqrt(Len*(Len-1)/2);

all_r = tril(all_r_R)+triu(all_r_S);
all_r = all_r - (all_r ==2);

all_r_abs = all_r .* (abs(all_r)>0.1);

abs_r_matrix =  xor(abs(all_r_S)>0.1, abs(all_r_R)>0.1);


all_r2_S = all_r_S.*all_r_S;
all_r2_R = all_r_R.*all_r_R;

%left CASE Right REFERENCE
all_r2 = tril(all_r2_R)+triu(all_r2_S);
all_r2 = all_r2 - (all_r2 ==2);

abs_r2_matrix = xor(abs(all_r2_S)>0.1, abs(all_r2_R)>0.1);

%If we define r2>0.1 as strong

sign_matrix = ((sign(all_r_S).*sign(all_r_R))==-1);
singDiffRate = sum(sum(sign_matrix))/(77*77);
plotWithSignSimple(all_r_R, all_r_S, zeros(77,77), 1);

plotPowerDist(all_r2);



end


function Tr = getTr(Y, r_S, r_R, maskMatrix, signMatrix)
if nargin == 3
    A2 = (2*Y'-1)*(2*Y-1);
    Tr = sum(sum((r_S - r_R).* A2))/2;
else
    r_S = abs(r_S).*signMatrix; %replace the sign
    r_S = r_S.*maskMatrix;
    r_R = r_R.*maskMatrix;
    A2 = (2*Y'-1)*(2*Y-1);
    Tr = sum(sum((r_S - r_R).* A2))/2;
end
end
