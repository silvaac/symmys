% this script familiarizes the user with some basic properties of the
% lognormal distribution, see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (1.94) and below

clear; close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs 
NumSimulations=10000; % scalar
mu=10;
sigma_2=.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate lognormal sample with above parameters
sigma=sqrt(sigma_2);
X=lognrnd(mu,sigma,NumSimulations,1);

% compute pdf and cdf using matlab functions
x=[min(X): (max(X)-min(X))/1000 : max(X)]; % set range
pdf =lognpdf(x,mu,sigma);
cdf =logncdf(x,mu,sigma);

% compute quantile using matlab functions
u=[.01 : .01 : .99]; % set range
quantile = logninv(u,mu,sigma);

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
h=plot(x,cdf,'b');
grid on
title('cdf')

figure
h=plot(u,quantile,'b');
grid on
title('quantile')