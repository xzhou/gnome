function [genoSeq] = hapSeq2GenoSeq(hapSeq, alleleMapping)
%hapSeq is fasta like haplotype sequence
[nHapSeq, nSnps] = size(hapSeq);
genoSeq = zeros(nHapSeq/2, nSnps);
% for i = 1:nHapSeq/2
%     for j = 1:nSnps
%         a1 = hapSeq(2*i,j);
%         a2 = hapSeq(2*i-1, j);
%         A = alleleMapping(j);
%         if a1 == A && a2 == A
%             genoSeq(i,j) = 0;
%         elseif a1 ~= A && a1 ~= A
%             genoSeq(i,j) = 2;
%         else
%             genoSeq(i,j) = 1;
%         end
%     end
% end

%%Here we begin the random select
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));


for i = 1: nHapSeq/2
    j = randi(nHapSeq, 1, 1);
    k = randi(nHapSeq, 1, 1);
    for l = 1:nSnps
        a1 = hapSeq(j, l);
        a2 = hapSeq(k, l);
        A = alleleMapping(l);
        if a1 == A && a2 == A
            genoSeq(i,l) = 0;
        elseif a1 ~= A && a1 ~= A
            genoSeq(i,l) = 2;
        else
            genoSeq(i,l) = 1;
        end
    end  
end
end