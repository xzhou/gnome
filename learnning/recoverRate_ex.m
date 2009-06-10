function [rate,ref]=recoverRate_ex(A,Aref,ratio)
if nargin == 2, ratio = 0.85; end
[row,col]=size(A);
[rowRef,colRef]=size(Aref);
occur(rowRef)=0;
used(row)=0;
ref(col) = 0;
for match=1.0:-0.05:ratio
    for i=1:row
        if used(i) == 1, continue; end
        for j=1:rowRef
            if occur(j) ~= 0, continue; end
            r = similarity(A(i,:),Aref(j,:));
            if r >= match
                a = A(i,:);
                aref = Aref(j,:);
                for n=1:length(A(i,:))
                    if a(n) ~= aref(n), 
                        ref(n) = ref(n)+1;
                    end
                end
                %fprintf(1,'rec:%s\n',a);
                occur(j) = r;
                used(i) = 1;
                break;
            end
        end
    end
end
rate=sum(occur)/length(occur);