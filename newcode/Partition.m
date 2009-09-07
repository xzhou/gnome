function [blocks] = partition(r, v)
global smap;
delta = 0.6;
smap = zeros(v,2);
S(r, v, delta);
blocks = smap;
end



function [value] = intraBlockSum(r, i, j, delta)
subMatrix = r(i:j, i:j) - delta;
value = sum(sum(subMatrix))/2;
end

function [value] = interBlockSum(r, i, j, k, delta)
leftSum = intraBlockSum(r, i, j-1, delta);
rightSum = intraBlockSum(r, j, k, delta);
total = intraBlockSum(r, i, k, delta);
value = total - leftSum -rightSum;
end


function [value, pos] = S(r, i, delta)
global smap;
if i ==0
    value = 0;
    pos = 0;
elseif i == 1
    value = 0;
    pos = 0;
elseif smap(i, 1) ~=0
    value = smap(i, 1);
    pos = smap(i, 2);
else
    allj = zeros(i-1, 1);
    for j=1:i-1
        [jv jp] = S(r, i-1, delta);
        intraj = intraBlockSum(r, j, i, delta);
        interj = interBlockSum(r, 1, j, i, delta);
        allj(j) = jv - interj + intraj;
    end
    [jv jpos] = max(allj(:,1));
    value = jv;
    pos = jpos;
    smap(i,:) = [value pos];
end
end
