% this script confirms that the analytical expression for the cdf of the gamma
% distribution is correct, see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% formula (1.111) 


clear; close all;  clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs 
nu=10;
sigma_square=1;

NumSimulations=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=nu/2;
b=2*sigma_square;

% generate sample using matlab functions
X=gamrnd(a,b,NumSimulations,1);

% compute pdf and cdf using matlab functions
x=[min(X): (max(X)-min(X))/1000 : max(X)]; % set range
cdf =gamcdf(x,a,b);

% lower regularized gamma function
F=gammainc(x/(2*sigma_square),nu/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
h=plot(x,cdf,'.');
hold on
h=plot(x,F,'r');
grid on