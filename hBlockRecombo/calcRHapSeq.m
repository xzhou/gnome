function [r] = calcRHapSeq(seq)
%%calculate r from 01 sequence, 1 as major, 0 as minor

%% commented for performance
% if ~isequal(sort(unique(seq)), [0 1])
%     e = MException('calcRHapSeq:error', 'seq must be 0 1 sequence');
%     throw(e);
% end

%% calculation
r = corrcoef(seq);
r(isnan(r)) = 0;
end