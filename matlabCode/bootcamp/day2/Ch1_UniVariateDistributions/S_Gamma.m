% this script familiarizes the user with the properties of the gamma
% distribution. It also shows how the chi-square distribution is a special
% case of the gamma distribution
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (1.108) and below

clear; close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs 
NumSimulations=10000; % scalar
nu=10;
sigma_square=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=nu/2;
b=2*sigma_square;

% generate sample using matlab functions
X=gamrnd(a,b,NumSimulations,1);

% compute pdf and cdf using matlab functions
x=[min(X): (max(X)-min(X))/1000 : max(X)]; % set range
pdf =gampdf(x,a,b);
cdf =gamcdf(x,a,b);
chiscdf =chi2cdf(x,nu);

% compute quantile using matlab functions
u=[.01 : .01 : .99]; % set range
quantile = gaminv(u,a,b);
chisquant = chi2inv(u,nu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
figure
NumBins=round(10*log(NumSimulations));
hist(X,NumBins)
[n,Ds]=hist(X,NumBins);
D=Ds(2)-Ds(1);
hold on
h=plot(x,pdf*NumSimulations*D,'r');
grid on
title('pdf/histogram')

figure
plot([1:NumSimulations],X,'.')
grid on
title('sample vs. observation time')

figure
h=plot(x,cdf,'g');
hold on
h=plot(x,chiscdf,'r');
grid on
legend('gamma','chi-square')
title('cdf')

figure
h=plot(u,quantile,'g');
hold on
h=plot(u,chisquant,'r');
grid on
legend('gamma','chi-square')
title('quantile')