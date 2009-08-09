function ID = getID(filename)
% get snp ID from list file
fid = fopen(filename,'r');
ID = textscan(fid,'%s');
ID = ID{1};