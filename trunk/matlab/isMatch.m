function ret=isMatch(a)

devi = a(1);
print = a(2);
a = a(3:end);
[row,col]=size(a);
y=rref3(a);
[r,c]=size(y);
min = 0;
max = 1;
for k=1:r
    if y(k,c-1) == -1
        tmpMin = -y(k,c);
        tmpMax = 1-y(k,c);
        if tmpMin > min
            min = tmpMin;
        end
        if tmpMax < max
            max = tmpMax;
        end    
    elseif y(k,c-1) == 1
        tmpMin = y(k,c)-1;
        tmpMax = y(k,c);
        if tmpMin > min
            min = tmpMin;
        end
        if tmpMax < max
            max = tmpMax;
        end
    end
end
if print ~= 0
    y
    max
    min
end
if (max-min) >= devi
    ret = 1;
else
    ret = 0;
end