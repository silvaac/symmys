% this script familiarizes the user with some basic properties of the
% Student t distribution, see "Risk and Asset Allocation"-Springer (2005), by A. Meucci
% Section 1.3.4

 clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
mu=0;
sigma_2=1;  
nu=1;  % 
NumSimulations=100000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute pdf and cdf using matlab functions
sigma=sqrt(sigma_2);
x=[mu-4*sigma : sigma/100 : mu+4*sigma]; % set range
pdf =1/sigma*tpdf((x-mu)/sigma,nu);
cdf =tcdf((x-mu)/sigma,nu);

% generate sample using matlab functions
X = mu+sigma*trnd(nu,NumSimulations,1);

% compute quantile using matlab functions
u=[.01 : .01 : .99]; % set range
quantile = mu+sigma*tinv(u,nu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
subplot(2,2,1)
NumBins=round(10*log(NumSimulations));
hist(X,NumBins)
[n,Ds]=hist(X,NumBins);
D=Ds(2)-Ds(1);
hold on
h=plot(x,pdf*NumSimulations*D,'r');
grid on
title('pdf/histogram')

subplot(2,2,2)
plot([1:NumSimulations],X,'.')
grid on
xlabel('realization number')
ylabel('realization value')

subplot(2,2,3)
h=plot(x,cdf,'b');
grid on
title('cdf')

subplot(2,2,4)
h=plot(u,quantile,'b');
grid on
title('quantile')