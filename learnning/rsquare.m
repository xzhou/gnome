function r=rsquare(x11,x12,x21,x22,precision)

if nargin == 4, precision = 1000; end

p1=x11+x12;
p2=x21+x22;
q1=x11+x21;
q2=x12+x22;
total = p1+p2;

div=p1*p2*q1*q2;
if div==0
    r = 0;
else
    r = (x11*x22-x12*x21)^2/div;
    %r = round(r * precision)/precision; 
end