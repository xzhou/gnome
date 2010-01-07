function [ seq ] = gInnerSeqLearning( targetSeq, refSeq, blocks, verbose )
%GINNERSEQLEARNING gInnerSeqLearning will return the sequence after inner
%block learning

  if nargin == 3
    verbose = false;
  end

  if verbose
    disp "staring genotype block learning"
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
  
  nRepeat = 10;
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
      fprintf(1, '\n************learning block %d ***************\n', i);
    end
    
    for k = 1:nRepeat
      [aResult] = gInnerBlockLearning(targetBlock, targetFreq, refBlock, refFreq);
      try
        result{i,k} = aResult;
      catch e
        fprintf(1, 'error\n')
      end
      if verbose
        fprintf(1, 'block = %d repeat = %d a = %f initSR = %f, finalSR = %f\n',i, k, aResult.fDistance, aResult.initSignRate, aResult.finalSignRate);
      end 
    end
  end
  
  %recover function will select the best quality result and reconstruct the
  %seqeunce
  
  seq = InnerBlockHelp.recoverCaseSeq(result);
  %save('innerseq.mat', 'seq');
  save;
  
end

