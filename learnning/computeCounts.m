function [x11,x12,x21,x22]=computeCounts(A,i,j)
[row,col]=size(A);
x11=0;
x12=0;
x21=0;
x22=0;
for k=1:row
    if A(k,i)=='0' && A(k,j)=='0'
        x11=x11+1;
    elseif A(k,i)=='0' && A(k,j)=='1'
        x12=x12+1;
    elseif A(k,i)=='1' && A(k,j)=='0'
        x21=x21+1;
    elseif A(k,i)=='1' && A(k,j)=='1'
        x22=x22+1;
    end
end