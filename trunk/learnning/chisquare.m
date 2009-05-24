function r=chisquare(x11,x12,x21,x22)
% compute chi-square
% x11, x12, x21, x22: four pairwise frequencies

p1=x11+x12;
p2=x21+x22;
q1=x11+x21;
q2=x12+x22;

div=p1*p2*q1*q2;
if div==0
    r = 0;
else
    r = ((x11*x22-x12*x21)^2)*(p1+p2)/div;
end