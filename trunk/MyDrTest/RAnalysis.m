cd '/home/xzhou/research_linux/gnome/workspace/data/dist100x77';
load;

numberOfR = zeros(0,3); %threshold | number of r values

step = 0.005;
k = 0;
for T = -1:step:0.95
    k = k + 1;
    tempR = R;
    index = logical((tempR > T) .* (tempR < (T + step)));
    
    tempR(index) = 1;
    tempR(~index) = 0;
    nr = sum(sum(tempR));
    
    nSignDiff = sum(sum(signDiff(index) == -1))/2;
    
    numberOfR(k, 1) = T;
    numberOfR(k, 2) = nr;
    numberOfR(k, 3) = nSignDiff;
end

absNumberOfR = zeros(0,3);
for T = 0:step:0.95
    k = k + 1;
    tempR = R;
    index = logical((abs(tempR) > T) .* (abs(tempR) < (T + step)));
    
    tempR(index) = 1;
    tempR(~index) = 0;
    nr = sum(sum(tempR));
    
    nSignDiff = sum(sum(signDiff(index) == -1))/2;
    
    absNumberOfR(k, 1) = T;
    absNumberOfR(k, 2) = nr;
    absNumberOfR(k, 3) = nSignDiff;
end

figure;
hold on;
plot(numberOfR(:,1), numberOfR(:,2),'color', 'red');
plot(numberOfR(:,1), numberOfR(:,3),'color', 'blue');
legend('r distribution', 'sign diff distribution');
saveas(gcf, 'r sign distribution.pdf');

figure;
hold on;
plot(absNumberOfR(:,1), absNumberOfR(:,2), 'color', 'green');
plot(absNumberOfR(:,1), absNumberOfR(:,3), 'color', 'black');
legend('r distribution', 'sign diff distribution');
saveas(gcf, 'abs r sign dist.pdf');

