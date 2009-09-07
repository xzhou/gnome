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

fre = zeros(77, 4);

for i = 1:Len
    [fre(i,1)  fre(i,2) fre(i,3) fre(i,4)] = getSingleAlleleFreq(int4S, i);
end

[nt1 fre1 nt2 fre2] = getSingleAlleleFreq(int4S, 64);




allele1 = int4S(end,:); %why use the last line as allele 0/1 definition
int2S = (int4S == repmat(allele1,nS,1)) + 0;
int2R = (int4R == repmat(allele1,nR,1)) + 0;

all_r_S = corrcoef(int2S);
all_r_S(isnan(all_r_S)) = 0;

all_r_R = corrcoef(int2R);
all_r_R(isnan(all_r_S)) = 0;

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

% 
 for i =1:sec
plotPowerDist(tempTrS(:,:,i));
 end

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

function [temp1 count1 temp2 count2] = getSingleAlleleFreq(int4S, col)
temp1 = int4S(1, col);
i = 1
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


