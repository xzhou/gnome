function i = randInt(min,max)

% randInt generates a random integer between min and max

i = fix((max-min+1)*rand(1)+min);