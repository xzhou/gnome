function genotypeBlock = innerBlockMutate(optimValues, problemData)
%This function difines the mutate function of the inner block learning

genotypeBlock = optimValues.x;
for i=1:floor(optimValues.temperature)+1
    [nrows ncols] = size(genotypeBlock);
    genotypeBlock = neighbor(genotypeBlock, nrows, ncols);
end
%=====================================================%


function genotypeBlock = neighbor(genotypeBlock, nrows, ncols)
row1 = randinteger(1,1,nrows)+1;
row2 = randinteger(1,1,nrows)+1;
while(genotypeBlock(row1,:)==genotypeBlock(row2,:))
    row1 = randinteger(1,1,nrows)+1;
    row2 = randinteger(1,1,nrows)+1;
end

genotypeBlock(row1,:) = genotypeBlock(row2,:);


%=====================================================%
function out = randinteger(m,n,range)
%RANDINTEGER generate integer random numbers (m-by-n) in range

len_range = size(range,1) * size(range,2);
% If the IRANGE is specified as a scalar.
if len_range < 2
    if range < 0
        range = [range+1, 0];
    elseif range > 0
        range = [0, range-1];
    else
        range = [0, 0];    % Special case of zero range.
    end
end
% Make sure RANGE is ordered properly.
range = sort(range);

% Calculate the range the distance for the random number generator.
distance = range(2) - range(1);
% Generate the random numbers.
r = floor(rand(m, n) * (distance+1));

% Offset the numbers to the specified value.
out = ones(m,n)*range(1);
out = out + r;