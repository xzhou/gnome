function [ ] = SIMDataAlleleMapingAnalysis( )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77';
    
    sampleFileName = 'SIM_100x77_smp.fasta';
    referenceFileName = 'SIM_100x77_ctl.fasta';
    
    seqS = fastaread(sampleFileName);
    seqR = fastaread(referenceFileName);
    
    mS = length(seqS);
    nS = length(seqS(1).Sequence);
    mR = length(seqR);
    nR = length(seqR(1).Sequence);
    
    int4S = zeros(mS, nS);
    int4R = zeros(mR, nR);
    
    for i = 1:mS
        int4S(i,:) = nt2int(seqS(i).Sequence) - 1;
    end
    
    for i = 1:mR
        int4R(i,:) = nt2int(seqR(i).Sequence) - 1;
    end
    
    smapping = getMajorAllele(int4S);
    rmapping = getMajorAllele(int4R);
    
    mixMapping = [smapping; rmapping];
    save('mixMapping', 'mixMapping');
    
    %NOTE: this is a random mapping, sepecial care must be taken to make it
    %correct
    mapping = int4S(end,:);
    int2S = (int4S == repmat(smapping, mS, 1)) + 0;
    int2R = (int4R == repmat(rmapping, mR, 1)) + 0;
    
    all_r_S = corrcoef(int2S);
    all_r_S(isnan(all_r_S)) = 0.0;
    all_r_R = corrcoef(int2R);
    all_r_R(isnan(all_r_R)) = 0.0;
    
    h = plotWithSignSimple(all_r_S, all_r_R, 0.1);
    saveas(h, 'mixMapping.pdf');
    
end


