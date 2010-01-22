function [ result ] = gInnerBlockLearning(caseBlock, caseFreq, refBlock, refFreq)
% using estimated r square value to learn learn genotype block
%% config
  expT = 1.0e-5;
  maxItr = 100;

%% learning
  targetGSeq = blockReconstruct(caseBlock, caseFreq);
  [targetR targetF targetC] = estimateR(targetGSeq);
  targetRs = targetR.*targetR;
  
  refGSeq = blockReconstruct(refBlock, refFreq);
  [refR refF refC] = estimateR(refGSeq);
  refRs = refR.*refR;
  
  currentFreq = refFreq;
  currentRs = refRs;
  currentF = refF;
  currentQuality = getQuality(targetRs, currentRs, targetF, currentF);
  currentR = refR;
  itr = 0;
  while itr < maxItr
    itr = itr + 1;
    [newSeq newFreq] = getNextSeq(refBlock, currentFreq);
    [newR newF newC] = estimateR(newSeq);
    newRs = newR.*newR;
    
    newQuality = getQuality(targetRs, newRs, targetF, newF);
    Qdiff = newQuality - currentQuality;
    
    % accepting probablity
    p = 1/(1+exp(Qdiff/expT));
    x = rand(1);
    if x < p
      currentFreq = newFreq;
      currentQuality = newQuality;
      currentRs = newRs;
      currentR = newR;
    end
  end
  
  result.finalFreq = currentFreq;
  result.finalQual = currentQuality;
  result.fianlR = currentR;
  result.refFreq = refFreq;
  result.refBlock = refBlock;
  result.initSignRate = InnerBlockHelp.calcSignRate(targetR, refR);
  result.finalSignRate = InnerBlockHelp.calcSignRate(targetR, currentR);
end