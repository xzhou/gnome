function haps=computeHaplotypes(A)
[row,col]=size(A);

haps = [];
for i=1:row
    [hRow,hCol]=size(haps);
    in = 0;
    for m=1:hRow
        if similarity(haps(m,:),A(i,:)) == 1
            in = 1;
            break;
        end
    end
    
    if in == 0
        haps = [haps;A(i,:)];
    end
end