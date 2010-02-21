function [] = printCoverRate(fid, blockCoverRate, str)
%function print type cover rate and sequence cover rate
fprintf(fid, '***%s***\n', str);
[tmp, nBlock] = size(blockCoverRate);
for i = 1:nBlock
    tr = blockCoverRate(i).typeCoverate;
    sr = blockCoverRate(i).seqCoverate;
    fprintf(fid, 'hyplotype block %d, type coverate = %f, seq coverrate = %f\n', i, tr, sr);
end
end