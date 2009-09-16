function test = RPower()

cd 'F:\IUBResearch\Projects\Bioinfor\data\dist100x77'

sample = 'SIM_100x77_smp.fasta';
reference = 'SIM_100x77_ctl.fasta';

%define the section length of r, for example sec = 20, 1/20 as the section
%length
sec = 20;

seqS = fastaread(sample);
seqR = fastaread(reference);

nS = length(seqS);
Len = length(seqS(1).Sequence);
nR = length(seqR);

int4S = zeros(nS, Len);
int4R = int4S;

for i = 1:length(seqS)
    int4S(i,:) = nt2int(seqS(i).Sequence) - 1;
    int4R(i,:) = nt2int(seqR(i).Sequence) - 1;
end

nt4S = zeros(nS, Len);
nt4R = nt4S;

nt4S = int2nt(int4S +1);
nt4R = int2nt(int4R + 1);

%For caculating the single allele frequency
freS = zeros(Len, 4);
freR= zeros(Len, 4);

for i = 1:Len
    [freS(i,1)  freS(i,2) freS(i,3) freS(i,4)] = getSingleAlleleFreq(int4S, i);
    [freR(i,1)  freR(i,2) freR(i,3) freR(i,4)] = getSingleAlleleFreq(int4R, i);
end

%For checking the allele frequence difference of the reference and case
%group
checkFre = zeros(Len, 1);

for i = 1:Len
    checkFre(i) = max(freS(i, 2), freS(i,4)) - max(freR(i,2), freR(i,4));
end
absfre = abs(checkFre);


%[nt1 fre1 nt2 fre2] = getSingleAlleleFreq(int4S, 64);


allele1 = int4S(end,:); %why use the last line as allele 0/1 definition
int2S = (int4S == repmat(allele1,nS,1)) + 0;
int2R = (int4R == repmat(allele1,nR,1)) + 0;

col1 = 17;
col2 = 23;

rS = 0;
rR = 0;
[rR rS] = getRValue(col1, col2, int4R, int4S);


all_r_S = corrcoef(int2S);
all_r_S(isnan(all_r_S)) = 0;

all_r_R = corrcoef(int2R);
all_r_R(isnan(all_r_S)) = 0;

%Get the number of SNPs pairs which are big on one side but small on other
%side and have different signs.

%diffRate is used for defining the difference between case and  reference
diffRate = 0.1;

% %Without filter
% diffMatrix = or(((abs(all_r_R)<diffRate).*(abs((all_r_S)>diffRate))), ((abs(all_r_S)<diffRate).*(abs(all_r_R)>diffRate)));
% %diffMatrix = ((abs(all_r_R)<diffRate).*(abs((all_r_S)>diffRate)));
% diffSignMatrix = or(((all_r_R<0).*(all_r_S>0)), ((all_r_S<0).*(all_r_R>0)));

%With filter(0.1)
% maskMatrix = abs(all_r_R)>0.1;
% filterR = maskMatrix.*all_r_R;
% filterS = maskMatrix.*all_r_S;
% diffMatrix = or(((abs(filterR)<diffRate).*(abs((filterS)>diffRate))), ((abs(filterS)<diffRate).*(abs(filterR)>diffRate)));
% diffSignMatrix = or(((filterR<0).*(filterS>0)), ((filterS<0).*(filterR>0)))

wholeMatrix = triu(all_r_R) + tril(all_r_S);
wholeMatrix = wholeMatrix - (wholeMatrix==2);
maskMatrix = abs(wholeMatrix)>0.1;
%diffMatrix = abs(all_r_R - all_r_S) > 0.25;
diffSignMatrix = or(((all_r_R<0).*(all_r_S>0)), ((all_r_S<0).*(all_r_R>0)));
%How to define the SNPs  pairs which are strong on one side but weak on the
%other side
diffSWMatrix = or(((abs(all_r_R)<diffRate).*(abs((all_r_S)>diffRate))), ((abs(all_r_S)<diffRate).*(abs(all_r_R)>diffRate)));
diffSMatrix = or(((abs(all_r_R)>diffRate).*(abs((all_r_S)>diffRate))), ((abs(all_r_S)>diffRate).*(abs(all_r_R)>diffRate)));
diffSWMatrix = diffSWMatrix.*maskMatrix;
diffSMatrix = diffSMatrix.*maskMatrix;

finalMatrix = wholeMatrix.*maskMatrix.*diffSignMatrix;
diffSignMatrix = diffSignMatrix.*maskMatrix;

signSWChangeRate =sum(sum(diffSignMatrix.*diffSWMatrix.*maskMatrix))/sum(sum(diffSWMatrix));
signSChanageRate = sum(sum(diffSignMatrix.*diffSMatrix.*maskMatrix))/sum(sum(diffSMatrix));



%For get temperory Tr by r distribution
StatS.Tr = zeros(nS,1);
StatR.Tr = zeros(nS,1);

averageS = zeros(50,1);
averageR = zeros(50,1);

tempTrS = zeros(Len, Len, sec);
tempTrR = zeros(Len, Len, sec);

for j = 1:sec
    %mask = or((j/sec-1/sec<=abs(all_r_S)).*(abs(all_r_S)<=j/sec), (j/sec-1/sec<=abs(all_r_R)).*(abs(all_r_R)<=j/sec));
    %both in reference and sample
    mask = (j/sec-1/sec<=abs(all_r_S)).*(abs(all_r_S)<=j/sec);
    %r only satisfy the requirement in Sample
    %mask = (j/sec-1/sec<=abs(all_r_R)).*(abs(all_r_R)<=j/sec);
    %r only satisfy the requirement in Reference
    maskS = all_r_S.*mask;
    maskR = all_r_R.*mask;

    %tempTrR is for calculating the power distribution according to R value
    for i = 1:nS
        StatS.Tr(i) = getTr(int2S(i,:), maskS, maskR);
        tempTrS(:,:,j) = tempTrS(:,:,j) + getTempTr(int2S(i,:), maskS, maskR);
        StatR.Tr(i) = getTr(int2R(i,:), maskS, maskR);
        tempTrR(:,:,j) = tempTrS(:,:,j) + getTempTr(int2R(i,:), maskS, maskR);
    end
    StatS.Tr = StatS.Tr/sqrt(Len*(Len-1)/2);  %??
    StatR.Tr = StatR.Tr/sqrt(Len*(Len-1)/2);
    averageS(j) = sum(StatS.Tr)/length(StatS.Tr);
    averageR(j) = sum(StatR.Tr)/length(StatR.Tr);
end

%Plotting without normalization
tempTrSAll = zeros(Len, Len, 1);

for i =1:sec
    tempTrSAll = tempTrS(:,:,i)+tempTrSAll;
    
end

plotPowerDist(tempTrSAll);

%py = 0:1/sec:1-1/sec;

%plot(averageR);
%plot(py,averageR,'k', py, averageS,'r');

% for i = 1:nS
%     StatS.Tr(i) = getTr(int2S(i,:), all_r_S, all_r_R);
%     StatR.Tr(i) = getTr(int2R(i,:), all_r_S, all_r_R);
% end
%
% StatS.Tr = StatS.Tr/sqrt(Len*(Len-1)/2);
% StatR.Tr = StatR.Tr/sqrt(Len*(Len-1)/2);
%
% plot(StatR.Tr, '.k');


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

function tempTr = getTempTr(Y, r_S, r_R)
A2 = (2*Y'-1)*(2*Y-1);
tempTr = (r_S - r_R).*A2;
%matrix for mapping the power of r by sections
end

function [h] = plotPowerDist(m, signMatrix)
[nrow ncol] = size(m);
h = figure();
hold on;
for i = 1:nrow
    for k = 1:ncol
        x_init = [k-1, k];
        y_init = [i-1, i];
        z_init = double([m(i,k) m(i,k); m(i,k) m(i,k)]);
        box = surface(x_init, y_init, z_init);
    end
end

maxm = max(max(abs(m)));

load('depColor', 'depColor');
set(gcf, 'ColorMap', depColor);
axis equal;
caxis([-maxm maxm]);
colorbar;

if nargin == 2
    signVector = MatrixToVec(signMatrix);
    nrow = length(signVector);
    z = ones(nrow, 1) * maxm + 0.1;

    scatter3(signVector(:,1), signVector(:,2), z,'marker', 'x');
end

h = gcf;
end

function  [rR rS] = getRValue(col1, col2, int4R, int4S)
[majorR1 minorR1] = getMajor(int4R, col1);
[majorR2  minorR2] = getMajor(int4R, col2);
[majorS1 minorS1] = getMajor(int4S, col1);
[majorS2  minorS2] = getMajor(int4S, col2);

C00R = 0;
C00S  = 0;
tempR = zeros(200, 2);
tempR(:,1) = int4R(:, col1);
tempR(:,2) = int4R(:, col2);
tempS = zeros(200, 2);
tempS(:,1) = int4S(:,col1);
tempS(:,2) = int4S(:, col2);
for i = 1:200
    if tempR(i, 1)  == majorR1 && tempR(i, 2) == majorR2
        C00R = C00R+1;
    end
end
for i = 1:200
    if tempS(i, 1) == majorS1 && tempS(i, 2) == majorS2
        C00S = C00S+1;
    end
end

freS = zeros(77, 4);
freR= zeros(77, 4);

for i = 1:77
    [freS(i,1)  freS(i,2) freS(i,3) freS(i,4)] = getSingleAlleleFreq(int4S, i);
    [freR(i,1)  freR(i,2) freR(i,3) freR(i,4)] = getSingleAlleleFreq(int4R, i);
end

rS = (C00S*200 - max(freS(col1, 2), freS(col1, 4))*max(freS(col2, 2), freS(col2, 4))*40000)/sqrt(freS(col1, 2)*freS(col1,4)*freS(col2, 2)*freS(col2,4)*40000*40000);
rR =  (C00R*200 - max(freR(col1, 2), freR(col1, 4))*max(freR(col2, 2), freR(col2, 4))*40000)/sqrt(freR(col1, 2)*freR(col1,4)*freR(col2, 2)*freR(col2,4)*40000*40000);
end

    


function [temp1 count1 temp2 count2] = getSingleAlleleFreq(int4S, col)
temp1 = int4S(1, col);
i = 1;
while(int4S(i,col) == temp1)
    i = i+1;
end
temp2 = int4S(i, col);



count1 = (sum(int4S(:,col)==temp1))/length(int4S);
count2 = (sum(int4S(:,col)~=temp1))/length(int4S);

switch temp1
    case 0
        temp1 = 'A'
    case 1
        temp1 = 'C'
    case 2
        temp1 = 'G'
    case 3
        temp1 = 'T'
end


switch temp2
    case 0
        temp2 = 'A'
    case 1
        temp2 = 'C'
    case 2
        temp2 = 'G'
    case 3
        temp2 = 'T'
end


% if count1>count2
%     major = count1;
%     minor = count2;
% else
%     major = count2;
%     minor = count1;
% end
end

function [major minor] = getMajor(int4S, col)
temp1 = int4S(1, col);
i = 1;
while(int4S(i,col) == temp1)
    i = i+1;
end
temp2 = int4S(i, col);

count1 = (sum(int4S(:,col)==temp1))/length(int4S);
count2 = (sum(int4S(:,col)~=temp1))/length(int4S);
if count1>count2
    major = temp1;
    minor = temp2;
else 
    major = temp2;
    minor = temp1;
end
end



