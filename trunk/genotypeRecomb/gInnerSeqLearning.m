function [ seq ] = gInnerSeqLearning(targetSeq, refSeq, config)
%GINNERSEQLEARNING gInnerSeqLearning will return the sequence after inner
%block learning

    verbose = config.verbose;
    nRepeat = config.nRepeat;
    blocks = config.blocks;
    
  if verbose
    disp 'staring genotype block learning'
  end
  
  [nBlock tmp] = size(blocks);
  [nS len] = size(targetSeq);
  
  targetBlockFreqInfo = cell(nBlock, 1);
  for i = 1:nBlock
    targetBlockFreqInfo{i,1} = gBlockFreq(targetSeq, blocks(i,:));
  end

  refBlockFreqInfo = cell(nBlock, 1);
  for i = 1:nBlock
    refBlockFreqInfo{i,1} = gBlockFreq(refSeq, blocks(i,:));
  end
  
  result = cell(nBlock, nRepeat);
  
  for i = 1:nBlock
    block = blocks(i,:);
    a = block(1,1);
    b = block(1,2);
    
    targetBlock = targetBlockFreqInfo{i,1}.uniqueBlock;
    targetFreq = targetBlockFreqInfo{i,1}.freq;
    blockTargetSeq = targetSeq(:, a:b);
    
    refBlock = refBlockFreqInfo{i,1}.uniqueBlock;
    refFreq = refBlockFreqInfo{i,1}.freq;
    blockRefSeq = refSeq(:,a:b);
    
    if verbose
      fprintf(config.logfid, '\n************learning block %d ***************\n', i);
    end
    %%TODO
    parfor k = 1:nRepeat
      [aResult] = gInnerBlockLearning(targetBlock, targetFreq, refBlock, refFreq, config);
      try
        result{i,k} = aResult;
      catch e
        fprintf(config.logfid, 'error\n')
      end
      if verbose
        fprintf(config.logfid, 'block = %d repeat = %d initQ = %f finalQ = %f initSR = %f, finalSR = %f\n',i, k, aResult.initQ, aResult.finalQual, aResult.initSignRate, aResult.finalSignRate);
      end 
    end
  end
  
  %recover function will select the best quality result and reconstruct the
  %seqeunce
  
  seq = InnerBlockHelp.recoverCaseSeq(result);
  %save('innerseq.mat', 'seq');
  save;
end

