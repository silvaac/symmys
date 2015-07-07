% this script shows how a generic distribution can be mapped into another arbitrary 
% distribution using the cdf-to-uniform + quantile-to-arbitrary trick
% see "Risk and Asset Allocation"-Springer (2005), by A. Meucci

close all; clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
NumSimulations=100000;
mu=-1;
sigma_2=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate normal and compute grade
sigma=sqrt(sigma_2);
X_a=normrnd(mu,sigma,NumSimulations,1);                        
U_a=normcdf(X_a,mu,sigma);

figure

NumBins=round(7*log(NumSimulations));
subplot(2,1,1)
hist(X_a,NumBins);
grid on 
xlabel('normal')

subplot(2,1,2)
hist(U_a,NumBins);
grid on 
xlabel('grade of normal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate uniform and anti-map into normal
U_b=rand(NumSimulations,1);
X_b=norminv(U_b,mu,sigma);

figure

subplot(2,1,2)
hist(U_b,NumBins);
grid on 
xlabel('uniform')

subplot(2,1,1)
hist(X_b,NumBins);
grid on 
xlabel('anti-grade of uniform')

