function rs=computeRSquaresForBlock(A,fun)

% A: a block of sequences
% fun: LD function (@dprime or @rsquare)
% rs: rsquares for all SNPs

if nargin == 1, fun=@rsquare; end  % fun in default is @rsquare

[row,col]=size(A);

for m=1:(col-1)
    for n=(m+1):col
        x11=0;
        x12=0;
        x21=0;
        x22=0;
        for i=1:row
            if A(i,m)=='0' && A(i,n)=='0'
                x11=x11+1;
            elseif A(i,m)=='0' && A(i,n)=='1'
                x12=x12+1;
            elseif A(i,m)=='1' && A(i,n)=='0'
                x21=x21+1;
            elseif A(i,m)=='1' && A(i,n)=='1'
                x22=x22+1;
            end
        end
        rs(m,n)=feval(fun,x11,x12,x21,x22);
        rs(n,m)=rs(m,n);
    end
end
[row,col]=size(rs);
for i=1:row
    rs(i,i)=sum(rs(:,i))/(row-1); % rs(i,i) has the averaged LD for ith SNP
end