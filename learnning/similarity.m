function rate=similarity(a,b)
%a sequence a
%b sequence b
%rate the fraction of SNPs with the same values in a and b
rate = 0;
if length(a) ~= length(b), 
    rate = 0;
else
    for i=1:length(a)
        if a(i) == b(i)
            rate = rate+1;
        end
    end
    rate = rate / length(a);
end