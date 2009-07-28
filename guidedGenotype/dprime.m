function r=dprime(x11,x12,x21,x22)

p1=x11+x12;
p2=x21+x22;
q1=x11+x21;
q2=x12+x22;
total = p1+p2;

D = x11*x22-x12*x21;
if D >=0
    Dmax = min(p1*q2,p2*q1);
else
    Dmax = max(-p1*q1,-p2*q2);
end
    
if Dmax == 0
    r = 0;
else
    r = D/Dmax;
end