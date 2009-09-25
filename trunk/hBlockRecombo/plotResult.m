function [test] = plotResult(index1, index2, preS, preR, postS, postR )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

figure;
hold on;

plot(index1, preS, '.r');
plot(index2, preR, '.g');

figure;
hold on;
plot(index1, postS, '.r');
plot(index2, postR, '.g');


end

