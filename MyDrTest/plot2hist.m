function power = plot2hist(D_M0,D_MA, textTitle,nbin, alpha)
% the function plot the histrogram of two datasets in one figure
% the hypothesis is that H0: T<=T0 vs. HA: T > T0;
if nargin <3, textTitle = {'histogram of two distribution';'H0 vs. HA'}; end
if nargin <4, nbin = 20; end
if nargin <5, alpha = 0.05; end % significance level
figure
min_D = min(min(D_M0), min(D_MA));
max_D = max(max(D_M0), max(D_MA));
binsize_D =  (max_D - min_D)/nbin;
bin = min_D:binsize_D:max_D;
n_M0 = hist(D_M0,bin);
n_MA = hist(D_MA,bin);
n_M0 = n_M0/sum(n_M0)/binsize_D;
n_MA = n_MA/sum(n_MA)/binsize_D;
plot(bin, n_M0, bin, n_MA);
legend({'H0: null' 'HA'});
title(textTitle);

% estimate the power;
D_M0 = sort(D_M0);
indx = floor(length(D_M0)*(1-alpha));
threshold = (D_M0(indx)+D_M0(indx+1))/2;
power = sum(D_MA>threshold)/length(D_MA);
end