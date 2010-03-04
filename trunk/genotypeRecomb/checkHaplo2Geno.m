function [] = checkHaplo2Geno(haploTypeSeq, wtccc1Conf);

    blocks = wtccc1Conf.blocks;
    [nBlock tmp] = size(blocks);

    haploTypeSeqFreqInfo = cell(nBlock, 1);
    
    parfor i = 1:nBlock
    haploTypeSeqFreqInfo{i,1} = getBlockFreq(haploTypeSeq, blocks(i,:));
    end

    [m n] = size(haploTypeSeq);
    
    for i = 1:nBlock
        [tempm, tempn] = size(haploTypeSeqFreqInfo{i,1});
        freq = 0;
        for j = 1:2:m
            k = 1;
            while sum(haploTypeSeq(j,blocks(i,1):blocks(i,2)) ~= haploTypeSeqFreqInfo{i, 1}(k, 1:end-1))~=0
                k = k+1;
            end
            if haploTypeSeqFreqInfo{i, 1}(k, end)>5
                freq1 = 1;
            else
                freq1 = 0;
            end
            
            k = 1;
            while sum(haploTypeSeq(j+1,blocks(i,1):blocks(i,2))~= haploTypeSeqFreqInfo{i, 1}(k, 1:end-1))~=0
                k = k+1;
            end
            if haploTypeSeqFreqInfo{i, 1}(k, end)>5
                freq2 = 1;
            else
                freq2 = 0;
            end
           
            freq = freq+freq1.*freq2;
        end
        perc = freq*2/m;
        fprintf(wtccc1Conf.logfid, 'Combination frequency between high Freq and hign Freq is %f for %d th block\n', perc, i);
    end
            
    
    for i = 1:nBlock
        [tempm, tempn] = size(haploTypeSeqFreqInfo{i,1});
        freq = 0;
        for j = 1:2:m           
            if sum(haploTypeSeq(j,blocks(i,1):blocks(i,2))~=haploTypeSeq(j+1,blocks(i,1):blocks(i,2)))==0
                freq = freq+1;
            end
        end
        perc = freq*2/m;
        fprintf(wtccc1Conf.logfid, 'Combination frequency between same Seq is %f for %d th block\n', perc, i);
    end
     

end