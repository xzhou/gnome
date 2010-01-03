function [result] = genotypeInterBlockRecomb(caseSeq, startingSeq, blocks, alleleMapping, config)
  %recombination of genotype sequence
  [m n] = size(blocks);
  
  if nargin == 4
    config.repeat = 100;
  end

  parfor nTrial = 1:config.repeat
    
  end
  
end
